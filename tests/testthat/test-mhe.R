test_that("unconstrained MHE equals the Kalman filter", {
  set.seed(4)
  sys <- ctrl_ss(matrix(c(1, 0.1, 0, 0.95), 2, byrow = TRUE), c(0, 0.1),
                 c(1, 0), dt = 0.1)
  Tn <- 40
  u <- sin(seq_len(Tn) / 5)
  x <- c(0.5, 0)
  y <- numeric(Tn)
  for (k in 1:Tn) {
    y[k] <- x[1] + rnorm(1, sd = 0.1)
    x <- drop(sys$A %*% x + sys$B * u[k]) + rnorm(2, sd = 0.03)
  }
  Qn <- diag(2) * 0.03^2
  kf <- ctrl_kalman_filter(sys, y, u, Qn = Qn, Rn = 0.01, x0 = c(0, 0),
                           P0 = diag(2))
  for (N in c(3L, 10L, 50L)) {
    mhe <- ctrl_mhe(sys, N, Qn = Qn, Rn = 0.01, x0 = c(0, 0), P0 = diag(2))
    est <- ctrl_mhe_run(mhe, y, u)
    expect_equal(unname(est$estimates), unname(kf$filtered), tolerance = 1e-6)
  }
})

test_that("constrained MHE respects state bounds where the Kalman filter does not", {
  set.seed(5)
  sys <- ctrl_ss(0.95, 1, 1, dt = 1)
  Tn <- 60
  x <- 0.05
  y <- numeric(Tn)
  xs <- numeric(Tn)
  for (k in 1:Tn) {
    xs[k] <- x
    y[k] <- x + rnorm(1, sd = 0.3)
    x <- max(0, 0.95 * x + rnorm(1, sd = 0.02))
  }
  kf <- ctrl_kalman_filter(sys, y, Qn = 4e-4, Rn = 0.09, x0 = 0, P0 = 1)
  mhe <- ctrl_mhe(sys, 10, Qn = 4e-4, Rn = 0.09, x0 = 0, P0 = 1, x_min = 0)
  est <- ctrl_mhe_run(mhe, y)
  expect_true(min(est$estimates) >= -1e-6)
  expect_true(min(kf$filtered) < 0)
  expect_lte(mean((est$estimates[, 1] - xs)^2), mean((kf$filtered[, 1] - xs)^2))
})

test_that("MHE refuses partially missing measurement rows", {
  sys <- ctrl_ss(diag(c(0.9, 0.8)), diag(2), diag(2), dt = 1)
  mhe <- ctrl_mhe(sys, 3, Qn = 0.01, Rn = 0.1)
  expect_error(ctrl_mhe_run(mhe, rbind(c(1, NA), c(1, 1))), "partially")
})

test_that("MHE handles singular arrival covariance and reports the full objective", {
  sys <- ctrl_ss(0, 1, 1, dt = 1)
  mhe <- ctrl_mhe(sys, 1, Qn = 1, Rn = 1, G = 0, P0 = 1)
  r <- ctrl_mhe_run(mhe, c(0, 0))
  expect_equal(drop(r$estimates), c(0, 0), tolerance = 1e-6)
  mhe1 <- ctrl_mhe(ctrl_ss(1, 0, 1, dt = 1), 1, Qn = 1, Rn = 1, x0 = 0, P0 = 1)
  r1 <- ctrl_mhe_run(mhe1, 2)
  expect_equal(r1$objective[1], 2, tolerance = 1e-4)
})
