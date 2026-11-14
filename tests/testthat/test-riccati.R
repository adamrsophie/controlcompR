test_that("scalar CARE matches sqrt(2) - 1", {
  r <- ctrl_care(-1, 1, 1, 1)
  expect_equal(drop(r$P), sqrt(2) - 1, tolerance = 1e-12)
  expect_equal(drop(r$K), sqrt(2) - 1, tolerance = 1e-12)
  expect_lt(r$residual, 1e-12)
})

test_that("scalar DARE matches the golden ratio", {
  r <- ctrl_dare(1, 1, 1, 1)
  P <- (1 + sqrt(5)) / 2
  expect_equal(drop(r$P), P, tolerance = 1e-12)
  expect_equal(drop(r$K), P / (1 + P), tolerance = 1e-12)
  expect_lt(Mod(r$poles), 1)
})
