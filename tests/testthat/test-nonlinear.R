pendulum <- function() {
  ctrl_nlsys(function(t, x, u, p) c(x[2], -p$g * sin(x[1]) - p$c * x[2] + u),
             output = function(t, x, u, p) x[1], nstates = 2, ninputs = 1,
             params = list(g = 9.81, c = 0.5))
}

test_that("linearization matches the analytic Jacobian", {
  sys <- pendulum()
  for (th in c(0, pi, 0.7)) {
    lin <- ctrl_linearize(sys, c(th, 0.3), 0.2)
    expect_equal(lin$A, matrix(c(0, 1, -9.81 * cos(th), -0.5), 2, byrow = TRUE),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(drop(lin$B), c(0, 1), tolerance = 1e-8)
    expect_equal(drop(lin$C), c(1, 0), tolerance = 1e-8)
  }
  up <- ctrl_linearize(sys, c(pi, 0))
  expect_false(ctrl_is_stable(up))
  expect_match(attr(up, "method"), "central")
})

test_that("equilibria and local stability", {
  sys <- pendulum()
  eq <- ctrl_equilibrium(sys, c(3, 0.1))
  expect_true(eq$converged)
  expect_equal(eq$x, c(pi, 0), tolerance = 1e-8)
  expect_equal(eq$stability, "unstable")
  down <- ctrl_equilibrium(sys, c(0.2, 0))
  expect_match(down$stability, "stable")
  # free input: hold the pendulum at 45 degrees
  eq2 <- ctrl_equilibrium(sys, c(pi / 4, 0), u0 = 0, free_inputs = 1,
                          fixed_states = 1)
  expect_equal(eq2$u, 9.81 * sin(pi / 4), tolerance = 1e-6)
})

test_that("nonlinear simulation matches linear simulation for a linear model", {
  A <- matrix(c(0, 1, -2, -0.4), 2, byrow = TRUE)
  lin <- ctrl_ss(A, c(0, 1), c(1, 0))
  nl <- ctrl_nlsys(function(t, x, u, p) drop(A %*% x) + c(0, u), nstates = 2,
                   ninputs = 1)
  t <- seq(0, 5, by = 0.1)
  u <- sin(t)
  a <- ctrl_nl_simulate(nl, c(1, 0), t, input = u)
  b <- ctrl_simulate(lin, input = u, times = t, x0 = c(1, 0))
  expect_equal(a$x[, 1], b$y[, 1, 1], tolerance = 1e-6)
  # discretized model reproduces the same samples
  nd <- ctrl_nl_discretize(nl, 0.1, method = "rk4", substeps = 10)
  c2 <- ctrl_nl_simulate(nd, c(1, 0), t, input = u)
  expect_equal(c2$x[, 1], b$y[, 1, 1], tolerance = 1e-6)
  # sampled feedback
  fb <- ctrl_nl_simulate(pendulum(), c(0.5, 0), seq(0, 10, by = 0.05),
                         controller = function(t, x) -5 * x[1] - 2 * x[2])
  expect_lt(abs(fb$x[201, 1]), 0.01)
})

test_that("EKF and UKF equal the Kalman filter on a linear model", {
  set.seed(9)
  A <- matrix(c(1, 0.1, 0, 0.9), 2, byrow = TRUE)
  lsys <- ctrl_ss(A, c(0, 0.1), c(1, 0), dt = 0.1)
  nl <- ctrl_nlsys(function(k, x, u, p) drop(A %*% x) + c(0, 0.1 * u),
                   output = function(k, x, u, p) x[1], nstates = 2, ninputs = 1,
                   dt = 0.1)
  y <- cumsum(rnorm(60, sd = 0.1))
  u <- rnorm(60)
  kf <- ctrl_kalman_filter(lsys, y, u, Qn = diag(c(1e-3, 1e-2)), Rn = 0.01)
  ekf <- ctrl_ekf(nl, y, u, Qn = diag(c(1e-3, 1e-2)), Rn = 0.01)
  ukf <- ctrl_ukf(nl, y, u, Qn = diag(c(1e-3, 1e-2)), Rn = 0.01)
  expect_equal(unname(ekf$filtered), unname(kf$filtered), tolerance = 1e-7)
  expect_equal(unname(ukf$filtered), unname(kf$filtered), tolerance = 1e-8)
  expect_equal(ukf$loglik, kf$loglik, tolerance = 1e-8)
})

test_that("EKF and UKF track a nonlinear pendulum", {
  set.seed(10)
  f <- function(k, x, u, p) c(x[1] + 0.05 * x[2], x[2] - 0.05 * 9.81 * sin(x[1]))
  h <- function(k, x, u, p) sin(x[1])
  sys <- ctrl_nlsys(f, h, nstates = 2, ninputs = 0, dt = 0.05)
  x <- c(1.2, 0)
  X <- matrix(0, 200, 2)
  y <- numeric(200)
  for (k in 1:200) {
    X[k, ] <- x
    y[k] <- sin(x[1]) + rnorm(1, sd = 0.05)
    x <- f(k, x, NULL, NULL) + rnorm(2, sd = 1e-3)
  }
  ekf <- ctrl_ekf(sys, y, Qn = 1e-6, Rn = 0.0025, x0 = c(1, 0), P0 = diag(c(0.1, 0.1)))
  ukf <- ctrl_ukf(sys, y, Qn = 1e-6, Rn = 0.0025, x0 = c(1, 0), P0 = diag(c(0.1, 0.1)))
  expect_lt(sqrt(mean((ekf$filtered[101:200, 1] - X[101:200, 1])^2)), 0.05)
  expect_lt(sqrt(mean((ukf$filtered[101:200, 1] - X[101:200, 1])^2)), 0.05)
})

test_that("NMPC agrees with linear MPC on a linear model and respects bounds", {
  A <- matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE)
  B <- c(0.005, 0.1)
  lsys <- ctrl_ss(A, B, diag(2), dt = 0.1)
  nl <- ctrl_nlsys(function(k, x, u, p) drop(A %*% x + B * u), nstates = 2,
                   ninputs = 1, dt = 0.1)
  mpc <- ctrl_mpc(lsys, 10, Q = diag(c(1, 0.1)), R = 0.1, u_min = -1, u_max = 1)
  nmpc <- ctrl_nmpc(nl, 10, Q = diag(c(1, 0.1)), R = 0.1, u_min = -1, u_max = 1)
  a <- ctrl_mpc_step(mpc, c(1, 0))
  b <- ctrl_nmpc_step(nmpc, c(1, 0))
  expect_true(b$ok)
  expect_equal(unname(b$u), unname(a$u), tolerance = 1e-3)
  sim <- ctrl_simulate_nmpc(nmpc, c(2, 0), 40)
  expect_true(all(abs(sim$u) <= 1 + 1e-12))
  expect_lt(abs(sim$x[41, 1]), 0.1)
})

test_that("EKF and UKF skip missing samples and use 0-based sample indices", {
  f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] - 0.1 * sin(x[1]))
  h <- function(k, x, u, p) x[1]
  sys <- ctrl_nlsys(f, h, nstates = 2, ninputs = 0, dt = 0.1)
  set.seed(1)
  x <- c(1, 0)
  y <- numeric(60)
  for (k in 1:60) {
    y[k] <- x[1] + rnorm(1, sd = 0.05)
    x <- f(k, x, NULL, NULL)
  }
  y[20:25] <- NA
  for (fn in list(ctrl_ekf, ctrl_ukf)) {
    r <- fn(sys, y, Qn = 1e-6, Rn = 0.0025, x0 = c(0.8, 0))
    expect_true(all(is.na(r$nis[20:25])))
    expect_true(all(is.finite(r$filtered)))
    expect_true(is.finite(r$loglik))
  }
  seen <- integer()
  f2 <- function(k, x, u, p) {
    seen <<- c(seen, k)
    x
  }
  sys2 <- ctrl_nlsys(f2, h, nstates = 2, ninputs = 0, dt = 0.1)
  ctrl_ukf(sys2, c(1, 1, 1), Qn = 1e-6, Rn = 1, x0 = c(1, 0))
  expect_equal(min(seen), 0L)
})

test_that("discrete outputs receive the sample index", {
  sys <- ctrl_nlsys(function(k, x, u, p) x, function(k, x, u, p) k,
                    nstates = 1, ninputs = 0, dt = 0.1)
  r <- ctrl_nl_simulate(sys, 0, times = c(0, 0.1, 0.2))
  expect_equal(drop(r$y), c(0, 1, 2))
  expect_error(ctrl_nl_simulate(sys, 0, times = c(0, 0.04, 0.08)), "integer multiples")
})
