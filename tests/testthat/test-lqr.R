test_that("continuous and discrete LQR match analytic references", {
  k <- ctrl_lqr(ctrl_ss(-1, 1, 1), 1, 1)
  expect_equal(drop(k$K), sqrt(2) - 1, tolerance = 1e-12)
  expect_equal(Re(k$poles), -sqrt(2), tolerance = 1e-12)
  kd <- ctrl_lqr(ctrl_ss(1, 1, 1, dt = 1), 1, 1)
  P <- (1 + sqrt(5)) / 2
  expect_equal(drop(kd$P), P, tolerance = 1e-12)
  expect_equal(drop(kd$K), P / (1 + P), tolerance = 1e-12)
  expect_equal(Re(kd$poles), 1 - P / (1 + P), tolerance = 1e-12)
})

test_that("closed loop is stable and has the LQR poles", {
  plant <- mass_spring()
  ctl <- ctrl_lqr(plant, diag(c(10, 1)), 1)
  cl <- ctrl_closed_loop(plant, ctl)
  expect_true(ctrl_is_stable(cl))
  expect_equal(sort_complex_test(ctrl_poles(cl)), sort_complex_test(ctl$poles),
               tolerance = 1e-10)
  expect_equal(nrow(cl$C), 2L)
  expect_error(ctrl_lqr(ctrl_ss(diag(c(1, -1)), c(0, 1), diag(2)), diag(2), 1),
               "stabilizable")
})

test_that("finite-horizon recursion converges to the DARE gain", {
  sys <- random_stable(4, 2, 1, seed = 5, dt = 0.1)
  inf <- ctrl_lqr(sys, diag(4), diag(2))
  fin <- ctrl_dlqr_finite(sys, diag(4), diag(2), horizon = 200)
  expect_equal(fin$K[[1]], unname(inf$K), tolerance = 1e-8)
  # with P_N = P_inf the finite horizon gain is exactly stationary
  fin2 <- ctrl_dlqr_finite(sys, diag(4), diag(2), horizon = 3, P_terminal = inf$P)
  expect_equal(fin2$K[[1]], unname(inf$K), tolerance = 1e-10)
})

test_that("reference gain gives unit DC gain", {
  plant <- mass_spring()
  ctl <- ctrl_lqr(plant, diag(2), 1)
  Nr <- ctrl_reference_gain(plant, ctl)
  cl <- ctrl_closed_loop(plant, ctl, include_control = FALSE)
  expect_equal(drop(ctrl_dcgain(cl) %*% Nr), 1, tolerance = 1e-12,
               ignore_attr = TRUE)
})
