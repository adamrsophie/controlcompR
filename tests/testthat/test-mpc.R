double_int <- function() {
  ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0)), 0.1)
}

test_that("unconstrained MPC with LQR terminal cost reproduces LQR", {
  plant <- double_int()
  Q <- diag(c(1, 0.5))
  R <- 0.2
  mpc <- ctrl_mpc(plant, horizon = 15, Q = Q, R = R, P = "lqr")
  K <- ctrl_lqr(plant, Q, R)$K
  for (x0 in list(c(1, 0), c(-0.3, 2))) {
    st <- ctrl_mpc_step(mpc, x0)
    expect_true(st$ok)
    expect_equal(st$u, drop(-K %*% x0), tolerance = 1e-6)
    expect_lt(st$dynamics_residual, 1e-7)
  }
})

test_that("input, rate and state constraints are respected", {
  plant <- double_int()
  mpc <- ctrl_mpc(plant, horizon = 20, Q = diag(c(10, 1)), R = 0.1,
                  u_min = -1, u_max = 1, du_min = -0.3, du_max = 0.3,
                  x_min = c(-Inf, -0.6), x_max = c(Inf, 0.6))
  sim <- ctrl_simulate_mpc(mpc, x0 = c(2, 0), steps = 80)
  expect_equal(sim$failures, 0)
  expect_true(all(abs(sim$u) <= 1 + 1e-6))
  expect_true(all(abs(diff(c(0, sim$u[, 1]))) <= 0.3 + 1e-6))
  expect_true(all(abs(sim$x[, 2]) <= 0.6 + 1e-5))
  expect_lt(abs(sim$x[81, 1]), 0.05)
  s <- summary(sim)
  expect_lt(s$max_input_violation, 1e-6)
})

test_that("control horizon holds the input", {
  plant <- double_int()
  mpc <- ctrl_mpc(plant, horizon = 10, control_horizon = 3, Q = diag(2), R = 1)
  st <- ctrl_mpc_step(mpc, c(1, 0))
  expect_equal(st$u_plan[4:10, 1], rep(st$u_plan[3, 1], 7), tolerance = 1e-7)
})

test_that("infeasible hard constraints are reported and soft constraints recover", {
  plant <- double_int()
  hard <- ctrl_mpc(plant, horizon = 5, Q = diag(2), R = 1, u_min = -0.1,
                   u_max = 0.1, x_min = c(-Inf, -0.05), x_max = c(Inf, 0.05))
  st <- ctrl_mpc_step(hard, c(0, 1))
  expect_false(st$ok)
  expect_true(all(is.na(st$u)))
  expect_error(ctrl_simulate_mpc(hard, c(0, 1), 3), "failed")
  sim <- ctrl_simulate_mpc(hard, c(0, 1), 3, on_failure = "hold")
  expect_true(all(sim$log$fallback == "hold"))
  soft <- ctrl_mpc(plant, horizon = 5, Q = diag(2), R = 1, u_min = -0.1,
                   u_max = 0.1, x_min = c(-Inf, -0.05), x_max = c(Inf, 0.05),
                   soft = "states")
  st2 <- ctrl_mpc_step(soft, c(0, 1))
  expect_true(st2$ok)
  expect_gt(max(st2$slack), 0)
  expect_lte(max(abs(st2$u)), 0.1 + 1e-6)
})

test_that("terminal equality constraint reaches the target", {
  plant <- double_int()
  mpc <- ctrl_mpc(plant, horizon = 30, Q = diag(2), R = 1, terminal = "equality")
  st <- ctrl_mpc_step(mpc, c(1, 0))
  expect_true(st$ok)
  expect_lt(max(abs(st$x_plan[31, ])), 1e-6)
})

test_that("objective reported equals direct evaluation", {
  plant <- double_int()
  Q <- diag(c(2, 1)); R <- 0.5
  mpc <- ctrl_mpc(plant, horizon = 6, Q = Q, R = R, u_min = -0.2, u_max = 0.2)
  st <- ctrl_mpc_step(mpc, c(1, 0.5))
  X <- st$x_plan; U <- st$u_plan
  J <- sum(vapply(1:6, function(k) sum(X[k, ] * (Q %*% X[k, ])) + R * U[k, 1]^2, 0)) +
    sum(X[7, ] * (Q %*% X[7, ]))
  expect_equal(st$objective, J, tolerance = 1e-6)
})

test_that("offset-free MPC removes steady-state error under mismatch", {
  model <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.1)
  truth <- ctrl_c2d(ctrl_tf(1.3, c(1, 1)), 0.1)
  mpc <- ctrl_mpc(model, 10, Q = 0, Qy = 1, R = 0.01, u_min = -5, u_max = 5)
  est <- ctrl_mpc_estimator(mpc, "output", Qn = 1e-4, Qd = 1e-3, Rn = 1e-4)
  sim <- ctrl_simulate_mpc(mpc, x0 = 0, steps = 200, plant = truth,
                           reference = list(y = 1), estimator = est,
                           disturbance = 0.2)
  expect_equal(sim$failures, 0)
  expect_equal(unname(sim$y[200, 1]), 1, tolerance = 1e-3)
})

test_that("MPC rejects invalid designs", {
  plant <- double_int()
  expect_error(ctrl_mpc(ctrl_ss(0, 1, 1), 5, 1, 1), "discrete")
  expect_error(ctrl_mpc(plant, 5, diag(2), 0), "positive definite")
  expect_error(ctrl_mpc(plant, 5, diag(c(1, -1)), 1), "semidefinite")
  expect_error(ctrl_mpc(plant, 5, diag(2), 1, u_min = 1, u_max = 0), "exceed")
})

test_that("simulation noise vectors are used per sample, not truncated", {
  plant <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.1)
  mpc <- ctrl_mpc(plant, 5, Q = 1, R = 0.1, u_min = -5, u_max = 5)
  est <- ctrl_mpc_estimator(mpc, "output", Qn = 1e-3, Qd = 1e-3, Rn = 1e-2)
  set.seed(3)
  nz <- rnorm(20, sd = 0.1)
  a <- ctrl_simulate_mpc(mpc, 0, 20, estimator = est, noise = nz,
                         reference = list(y = 1))
  b <- ctrl_simulate_mpc(mpc, 0, 20, estimator = est, noise = nz[1],
                         reference = list(y = 1))
  expect_false(isTRUE(all.equal(a$u, b$u)))
  expect_error(ctrl_simulate_mpc(mpc, 0, 20, noise = c(1, 2)), "length")
})

test_that("LQR terminal weight includes the output weight", {
  plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                            c(1, 0)), 0.1)
  m <- ctrl_mpc(plant, 10, Q = 0, Qy = 2, R = 0.2, P = "lqr")
  K <- ctrl_lqr(plant, diag(c(2, 0)), 0.2)$K
  expect_equal(unname(ctrl_mpc_step(m, c(1, 0))$u), unname(drop(-K %*% c(1, 0))),
               tolerance = 1e-6)
})

test_that("output bounds cover the applied input through feedthrough", {
  plant <- ctrl_ss(0, 0, 0, D = 1, dt = 1)
  m <- ctrl_mpc(plant, 2, Q = 0, Qy = 1, R = 1, y_max = 1)
  st <- ctrl_mpc_step(m, 0, reference = list(y = 10))
  expect_lte(unname(st$u), 1 + 1e-4)
  st2 <- ctrl_mpc_step(m, 0, reference = list(y = 10), output_offset = 0.5)
  expect_lte(unname(st2$u), 0.5 + 1e-4)
})

test_that("slack weights must be nonnegative", {
  plant <- ctrl_ss(1, 1, 1, dt = 1)
  expect_error(ctrl_mpc(plant, 3, 1, 1, x_max = 1, soft = "states",
                        soft_weight = c(-1, 0)), "nonnegative")
})
