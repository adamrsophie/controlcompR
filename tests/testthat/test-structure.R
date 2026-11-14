
test_that("dc gain ignores hidden poles at the origin; place takes poles positionally", {
  expect_equal(unname(ctrl_dcgain(ctrl_ss(diag(c(0, -1)), c(0, 1), c(1, 1)))), matrix(1))
  expect_error(ctrl_dcgain(ctrl_tf(1, c(1, 0))), "undefined")
  plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
  K <- ctrl_place(plant, c(-4, -5))
  expect_equal(sort(Re(eigen(plant$A - plant$B %*% K)$values)), c(-5, -4))
})

test_that("controllability survives overflowing matrix powers", {
  res <- ctrl_controllability(1000 * diag(105), diag(105))
  expect_true(res$controllable)
})
