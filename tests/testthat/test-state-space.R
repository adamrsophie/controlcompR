test_that("ctrl_ss validates dimensions and time domain", {
  sys <- mass_spring()
  expect_s3_class(sys, "ctrl_ss")
  expect_equal(unname(ctrl_dims(sys)), c(2L, 1L, 1L))
  expect_true(ctrl_is_continuous(sys))
  expect_equal(sys$D, matrix(0, 1, 1))
  expect_error(ctrl_ss(diag(2), c(1, 1, 1), c(1, 0)), "B")
  expect_error(ctrl_ss(diag(2), c(0, 1), c(1, 0, 0)), "`C`")
  expect_error(ctrl_ss(diag(2), c(0, 1), c(1, 0), D = matrix(0, 2, 2)), "1 x 1")
  expect_error(ctrl_ss(matrix(c(1, NA, 0, 1), 2), c(0, 1), c(1, 0)), "non-finite")
  expect_error(ctrl_ss(-1, 1, 1, dt = -0.1), "dt")
  expect_error(ctrl_ss(-1, 1, 1, dt = NA), "dt")
  expect_error(ctrl_ss(matrix(1:6, 2), 1, 1), "square")
  d <- ctrl_ss(0.5, 1, 1, dt = 0.1)
  expect_true(ctrl_is_discrete(d))
  expect_equal(ctrl_sample_time(d), 0.1)
})

test_that("signal names, subsetting and static gains", {
  sys <- ctrl_ss(-diag(2), diag(2), diag(2), input_names = c("a", "b"),
                 output_names = c("p", "q"))
  sub <- sys["q", "a"]
  expect_equal(sub$B, matrix(c(1, 0), 2))
  expect_equal(sub$output_names, "q")
  expect_error(ctrl_ss(-1, 1, 1, input_names = c("a", "b")), "length 1")
  expect_error(ctrl_ss(-diag(2), diag(2), diag(2), input_names = c("a", "a")),
               "unique")
  g <- ctrl_gain(matrix(1:4, 2))
  expect_equal(unname(ctrl_dims(g)), c(0L, 2L, 2L))
  expect_equal(ctrl_dcgain(g), matrix(1:4, 2), ignore_attr = TRUE)
})

test_that("similarity transforms preserve the transfer function", {
  sys <- mass_spring()
  T <- matrix(c(2, 1, 0, 1), 2)
  st <- ctrl_transform(sys, T)
  w <- c(0.1, 1, 10)
  expect_equal(ctrl_freqresp(st, w)$response, ctrl_freqresp(sys, w)$response,
               tolerance = 1e-12)
  expect_output(print(sys), "State-space model")
  expect_output(print(summary(sys)), "asymptotically stable")
})
