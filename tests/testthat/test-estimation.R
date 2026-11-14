test_that("ctrl_place assigns SISO and MIMO poles", {
  sys <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
  K <- ctrl_place(sys, poles = c(-4, -5))
  ev <- sort(Re(eigen(sys$A - sys$B %*% K)$values))
  expect_equal(ev, c(-5, -4), tolerance = 1e-8)
  set.seed(3)
  A <- matrix(rnorm(16), 4)
  B <- matrix(rnorm(8), 4)
  p <- c(-1, -2, -3 + 1i, -3 - 1i)
  K <- ctrl_place(A, B, p)
  ach <- eigen(A - B %*% K)$values
  expect_lt(max(vapply(p, function(z) min(Mod(ach - z)), 0)), 1e-7)
  expect_true(is.finite(attr(K, "cond_X")))
  expect_error(ctrl_place(A, B, c(-1, -1, -1, -2)), "repeated")
  expect_error(ctrl_place(A, B, c(-1, -2, -3 + 1i, -4)), "conjugate")
})

test_that("ctrl_place refuses uncontrollable pairs", {
  A <- diag(c(1, 2))
  B <- matrix(c(1, 0))
  expect_error(ctrl_place(A, B, c(-1, -2)), "not controllable")
})

test_that("observer poles are placed", {
  sys <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
  obs <- ctrl_observer(sys, c(-10, -11))
  ev <- sort(Re(eigen(sys$A - obs$L %*% sys$C)$values))
  expect_equal(ev, c(-11, -10), tolerance = 1e-8)
})

test_that("steady-state Kalman filter matches the limit of the time-varying filter", {
  sys <- ctrl_ss(matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE), c(0.005, 0.1),
                 c(1, 0), dt = 0.1)
  Qn <- diag(c(1e-4, 1e-3))
  kf <- ctrl_kalman(sys, Qn = Qn, Rn = 0.01)
  P <- diag(2)
  for (k in 1:2000) {
    S <- drop(sys$C %*% P %*% t(sys$C)) + 0.01
    M <- P %*% t(sys$C) / S
    P <- sys$A %*% (P - M %*% sys$C %*% P) %*% t(sys$A) + Qn
  }
  expect_equal(kf$P, P, tolerance = 1e-8)
  expect_true(all(Mod(kf$poles) < 1))
  expect_equal(unname(kf$L), unname(sys$A %*% kf$M), tolerance = 1e-10)
  expect_error(ctrl_kalman(sys, Qn, Rn = 0), "positive definite")
})

test_that("continuous Kalman-Bucy filter is stabilizing", {
  sys <- ctrl_ss(matrix(c(0, 1, -1, -0.2), 2, byrow = TRUE), c(0, 1), c(1, 0))
  kf <- ctrl_kalman(sys, Qn = diag(2), Rn = 0.1)
  expect_true(all(Re(kf$poles) < 0))
  res <- sys$A %*% kf$P + kf$P %*% t(sys$A) + diag(2) -
    kf$P %*% t(sys$C) %*% sys$C %*% kf$P / 0.1
  expect_lt(max(abs(res)), 1e-9)
})

test_that("time-varying filter converges to steady-state gain and smoother reduces covariance", {
  set.seed(1)
  sys <- ctrl_ss(matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE), c(0.005, 0.1),
                 c(1, 0), dt = 0.1)
  Qn <- diag(c(1e-6, 1e-4))
  x <- c(0, 1)
  y <- numeric(300)
  for (k in 1:300) {
    y[k] <- x[1] + rnorm(1, sd = 0.1)
    x <- drop(sys$A %*% x) + rnorm(2, sd = sqrt(diag(Qn)))
  }
  kf <- ctrl_kalman_filter(sys, y, Qn = Qn, Rn = 0.01)
  ss <- ctrl_kalman(sys, Qn = Qn, Rn = 0.01)
  expect_equal(kf$gain[, , 300], drop(ss$M), tolerance = 1e-5)
  expect_true(is.finite(kf$loglik))
  sm <- ctrl_kalman_smoother(kf)
  expect_true(all(diag(sm$P_smoothed[, , 150]) <= diag(kf$P_filtered[, , 150]) + 1e-12))
  y2 <- y
  y2[10] <- NA
  kf2 <- ctrl_kalman_filter(sys, y2, Qn = Qn, Rn = 0.01)
  expect_true(is.na(kf2$nis[10]))
})

test_that("LQG closed-loop poles are the union of regulator and estimator poles", {
  plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                            c(1, 0)), 0.1)
  for (form in c("current", "delayed")) {
    lqg <- ctrl_lqg(plant, Q = diag(c(1, 0.1)), R = 0.1, Qn = diag(c(1e-4, 1e-3)),
                    Rn = 1e-2, form = form)
    cl <- ctrl_closed_loop(plant, lqg)
    got <- sort_complex_test(ctrl_poles(cl))
    want <- sort_complex_test(c(lqg$lqr$poles, lqg$kalman$poles))
    expect_equal(got, want, tolerance = 1e-6)
  }
  csys <- ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0))
  lqg <- ctrl_lqg(csys, Q = diag(2), R = 1, Qn = diag(2), Rn = 1)
  got <- sort_complex_test(ctrl_poles(ctrl_closed_loop(csys, lqg)))
  want <- sort_complex_test(c(lqg$lqr$poles, lqg$kalman$poles))
  expect_equal(got, want, tolerance = 1e-6)
})

test_that("stochastic LQG simulation is reproducible", {
  plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                            c(1, 0)), 0.1)
  lqg <- ctrl_lqg(plant, Q = diag(c(1, 0.1)), R = 0.1, Qn = diag(c(1e-4, 1e-3)),
                  Rn = 1e-2)
  a <- ctrl_simulate_lqg(plant, lqg, 50, c(1, 0), replications = 3, seed = 7)
  b <- ctrl_simulate_lqg(plant, lqg, 50, c(1, 0), replications = 3, seed = 7)
  expect_identical(a$x, b$x)
  expect_lt(abs(a$x[51, 1, 1]), 0.5)
})

test_that("LQG simulation draws correlated noise and costs the cross weight", {
  plant <- ctrl_ss(0.9, 1, 1, dt = 1)
  lqg <- ctrl_lqg(plant, Q = 1, R = 1, Qn = 1, Rn = 1, Nn = 0.5, form = "delayed")
  sim <- ctrl_simulate_lqg(plant, lqg, steps = 4000, x0 = 0, seed = 3)
  x <- sim$x[, 1, 1]
  y <- sim$y[, 1, 1]
  w <- x[-1] - 0.9 * x[-length(x)] - sim$u[, 1, 1]
  v <- y - x[-length(x)]
  expect_equal(cov(w, v), 0.5, tolerance = 0.1)
  expect_error(ctrl_simulate_lqg(plant, lqg, steps = 3, x0 = 0, Qn = -1), "semidefinite")
  expect_error(ctrl_kalman(ctrl_ss(-1, 1, 1), Qn = 0, Rn = 1, Nn = 1, G = 0),
               "semidefinite")
})

test_that("finite-horizon LQR rejects indefinite weights", {
  expect_error(ctrl_dlqr_finite(ctrl_ss(1, 1, 1, dt = 1), 1, 1, 1, P_terminal = -2),
               "semidefinite")
})

test_that("RTS smoother accepts singular covariances", {
  kf <- ctrl_kalman_filter(ctrl_ss(1, 1, 1, dt = 1), y = c(0, 0), Qn = 0, Rn = 1, P0 = 0)
  sm <- ctrl_kalman_smoother(kf)
  expect_equal(drop(sm$smoothed), c(0, 0))
})
