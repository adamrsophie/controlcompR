test_that("zero-order hold is exact for scalar, singular and unstable systems", {
  for (Ts in c(0.01, 0.1, 1)) {
    d <- ctrl_c2d(ctrl_ss(-1, 1, 1), Ts)
    expect_equal(drop(d$A), exp(-Ts), tolerance = 1e-14)
    expect_equal(drop(d$B), 1 - exp(-Ts), tolerance = 1e-14)
    u <- ctrl_c2d(ctrl_ss(2, 1, 1), Ts)
    expect_equal(drop(u$A), exp(2 * Ts), tolerance = 1e-13)
    expect_equal(drop(u$B), (exp(2 * Ts) - 1) / 2, tolerance = 1e-13)
    dbl <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0)), Ts)
    expect_equal(dbl$A, matrix(c(1, 0, Ts, 1), 2), tolerance = 1e-14)
    expect_equal(drop(dbl$B), c(Ts^2 / 2, Ts), tolerance = 1e-14)
  }
  expect_error(ctrl_c2d(ctrl_ss(-1, 1, 1, dt = 1), 0.1), "already")
})

test_that("first-order hold, Tustin and inverse maps", {
  # integrator with FOH: x(k+1) = x(k) + T (u(k) + u(k+1)) / 2
  f <- ctrl_c2d(ctrl_ss(0, 1, 1), 0.5, "foh")
  expect_equal(drop(f$A), 1)
  r <- ctrl_simulate(f, input = 0:4, times = (0:4) * 0.5)
  expect_equal(r$y[, 1, 1], ((0:4) * 0.5)^2, tolerance = 1e-12)
  tu <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.1, "tustin")
  # bilinear: (T/2)(z+1) / ((1 + T/2) z - (1 - T/2))
  num <- c(0.05, 0.05) / 1.05
  den <- c(1, -0.95 / 1.05)
  expect_equal(tu$num[[1]], num, tolerance = 1e-12)
  expect_equal(tu$den[[1]], den, tolerance = 1e-12)
  sys <- random_stable(3, 2, 2, seed = 4)
  w <- c(0.3, 1, 3)
  for (m in c("zoh", "tustin")) {
    back <- ctrl_d2c(ctrl_c2d(sys, 0.05, m), m)
    expect_equal(ctrl_freqresp(back, w)$response, ctrl_freqresp(sys, w)$response,
                 tolerance = 1e-8)
  }
  pw <- ctrl_c2d(ctrl_tf(1, c(1, 0.1, 4)), 0.3, "tustin", prewarp = 2)
  gc <- ctrl_freqresp(ctrl_tf(1, c(1, 0.1, 4)), 2)$response
  gd <- ctrl_freqresp(pw, 2)$response
  expect_equal(gd, gc, tolerance = 1e-10)
})

test_that("matched pole-zero mapping preserves DC gain and maps poles", {
  G <- ctrl_tf(c(1, 2), c(1, 4, 3))
  Gd <- ctrl_c2d(G, 0.1, "matched")
  expect_equal(sort(Re(ctrl_poles(Gd))), sort(exp(c(-1, -3) * 0.1)), tolerance = 1e-12)
  expect_equal(drop(ctrl_dcgain(Gd)), drop(ctrl_dcgain(G)), tolerance = 1e-12)
})
