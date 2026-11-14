test_that("poles and stability classification", {
  expect_equal(sort(Re(ctrl_poles(mass_spring()))), c(-2, -1), tolerance = 1e-12)
  expect_equal(ctrl_stability(mass_spring())$classification, "asymptotically stable")
  dbl <- ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0))
  expect_equal(ctrl_stability(dbl)$classification, "unstable")
  expect_true(ctrl_stability(dbl)$defective_boundary)
  osc <- ctrl_ss(matrix(c(0, 1, -1, 0), 2, byrow = TRUE), c(0, 1), c(1, 0))
  expect_equal(ctrl_stability(osc)$classification, "marginally stable")
  # repeated semisimple eigenvalue at 0 is marginally stable
  expect_equal(ctrl_stability(ctrl_ss(matrix(0, 2, 2), diag(2), diag(2)))$classification,
               "marginally stable")
  expect_equal(ctrl_stability(ctrl_ss(1, 1, 1))$classification, "unstable")
  expect_equal(ctrl_stability(ctrl_ss(0.5, 1, 1, dt = 1))$classification,
               "asymptotically stable")
  expect_equal(ctrl_stability(ctrl_ss(-1, 1, 1, dt = 1))$classification,
               "marginally stable")
  expect_equal(ctrl_stability(ctrl_ss(1.01, 1, 1, dt = 1))$classification, "unstable")
  expect_equal(ctrl_stability(ctrl_tf(1, c(1, 0, 0)))$classification, "unstable")
})

test_that("damping, DC gain and zeros", {
  d <- ctrl_damp(ctrl_tf(1, c(1, 0.4, 4)))
  expect_equal(d$frequency, c(2, 2), tolerance = 1e-12)
  expect_equal(d$damping, c(0.1, 0.1), tolerance = 1e-12)
  expect_equal(drop(ctrl_dcgain(ctrl_tf(2, c(1, 3, 2)))), 1, tolerance = 1e-12)
  expect_error(ctrl_dcgain(ctrl_tf(1, c(1, 0))), "pole")
  dd <- ctrl_c2d(ctrl_tf(1, c(1, 0.4, 4)), 0.05)
  expect_equal(ctrl_damp(dd)$damping, c(0.1, 0.1), tolerance = 1e-8)
  z <- ctrl_zeros(ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(3, 1)))
  expect_equal(Re(z), -3, tolerance = 1e-10)
  expect_equal(Re(ctrl_zeros(ctrl_tf(c(1, 2), c(1, 3, 2)))), -2, tolerance = 1e-12)
})
