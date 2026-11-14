test_that("free response of x' = -x is exp(-t) on regular and irregular grids", {
  sys <- ctrl_ss(-1, 1, 1)
  t1 <- seq(0, 5, by = 0.1)
  expect_equal(ctrl_simulate(sys, times = t1, x0 = 1)$x[, 1, 1], exp(-t1),
               tolerance = 1e-13)
  t2 <- c(0, 0.01, 0.3, 0.31, 1, 2.5, 4)
  expect_equal(ctrl_initial(sys, 1, t2)$y[, 1, 1], exp(-t2), tolerance = 1e-13)
})

test_that("step and impulse responses match analytic solutions", {
  sys <- ctrl_ss(-1, 1, 1)
  t <- seq(0, 6, by = 0.25)
  expect_equal(ctrl_step(sys, t)$y[, 1, 1], 1 - exp(-t), tolerance = 1e-13)
  expect_equal(ctrl_impulse(sys, t)$y[, 1, 1], exp(-t), tolerance = 1e-13)
  d <- ctrl_ss(0.5, 1, 1, dt = 1)
  expect_equal(ctrl_impulse(d, 0:5)$y[, 1, 1], c(0, 0.5^(0:4)))
  expect_equal(ctrl_step(d, 0:4)$y[, 1, 1], c(0, 1, 1.5, 1.75, 1.875))
  expect_warning(ctrl_impulse(ctrl_ss(-1, 1, 1, D = 1), t), "Dirac")
  mimo <- ctrl_step(random_stable(3, 2, 2, seed = 9))
  expect_equal(dim(mimo$y)[3], 2L)
})

test_that("forced response with ramp input is exact with FOH", {
  sys <- ctrl_ss(0, 1, 1)  # integrator
  t <- seq(0, 2, by = 0.2)
  r <- ctrl_simulate(sys, input = t, times = t, hold = "foh")
  expect_equal(r$y[, 1, 1], t^2 / 2, tolerance = 1e-13)
  r0 <- ctrl_simulate(sys, input = t, times = t, hold = "zoh")
  expect_equal(r0$y[, 1, 1], c(0, cumsum(t[-length(t)] * 0.2)), tolerance = 1e-13)
  f <- ctrl_simulate(sys, input = function(tt) 1, times = t)
  expect_equal(f$y[, 1, 1], t, tolerance = 1e-13)
})

test_that("discrete simulation checks the time grid and input shape", {
  d <- ctrl_ss(0.5, 1, 1, dt = 0.1)
  expect_error(ctrl_simulate(d, input = rep(1, 3), times = c(0, 0.1, 0.25)),
               "sample instants")
  expect_error(ctrl_simulate(ctrl_ss(-1, 1, 1), input = 1:3, times = 1:4), "one row")
  df <- as.data.frame(ctrl_step(d, (0:3) * 0.1))
  expect_true(all(c("time", "case", "kind", "signal", "value") %in% names(df)))
})

test_that("default step horizon covers slow poles of stiff models", {
  G <- ctrl_tf(1, c(0.1, 100.001, 1))
  r <- ctrl_step(G)
  expect_gt(max(r$time), 500)
  expect_false(is.na(ctrl_step_info(G)$settling_time))
})
