test_that("normalization, printing and properness", {
  G <- ctrl_tf(c(0, 2, 4), c(2, 6, 4))
  expect_equal(G$num[[1]], c(1, 2))
  expect_equal(G$den[[1]], c(1, 3, 2))
  expect_output(print(G), "s\\^2")
  expect_false(ctrl_is_proper(ctrl_tf(c(1, 0, 0), c(1, 1))))
  expect_error(ctrl_tf2ss(ctrl_tf("s")), "improper")
  expect_error(ctrl_tf(1, 0), "zero")
  expect_error(ctrl_tf("z"), "dt")
})

test_that("tf and ss round trips preserve frequency responses", {
  G <- ctrl_tf(c(2, 3, 1), c(1, 4, 6, 4))
  S <- ctrl_tf2ss(G)
  w <- c(0.01, 0.5, 2, 20)
  expect_equal(ctrl_freqresp(S, w)$response, ctrl_freqresp(G, w)$response,
               tolerance = 1e-12)
  back <- ctrl_ss2tf(S)
  expect_equal(back$num[[1]], G$num[[1]], tolerance = 1e-10)
  expect_equal(back$den[[1]], G$den[[1]], tolerance = 1e-10)
  # biproper
  Gb <- ctrl_tf(c(3, 1), c(1, 2))
  expect_equal(ctrl_tf2ss(Gb)$D, matrix(3))
  mimo <- random_stable(3, 2, 2, seed = 7)
  tfm <- ctrl_ss2tf(mimo)
  expect_equal(ctrl_freqresp(tfm, w)$response, ctrl_freqresp(mimo, w)$response,
               tolerance = 1e-8)
  ssm <- ctrl_tf2ss(tfm)
  expect_equal(ctrl_freqresp(ssm, w)$response, ctrl_freqresp(mimo, w)$response,
               tolerance = 1e-8)
  expect_lte(ctrl_dims(ssm)[["states"]], 6L)
})

test_that("zpk conversions and minimal realization", {
  H <- ctrl_zpk(-1, c(-2, -3 + 1i, -3 - 1i), 4)
  Ht <- ctrl_zpk2tf(H)
  expect_equal(Ht$num[[1]], c(4, 4))
  back <- ctrl_tf2zpk(Ht)
  expect_equal(sort(Re(back$poles)), c(-3, -3, -2), tolerance = 1e-10)
  expect_error(ctrl_zpk(poles = c(-1 + 1i)), "conjugate")
  cancel <- ctrl_tf_minreal(ctrl_tf(c(1, 1), c(1, 3, 2)))
  expect_equal(cancel$den[[1]], c(1, 2), tolerance = 1e-10)
  red <- ctrl_minreal(ctrl_ss(diag(c(-1, -1)), c(1, 1), c(1, 1)))
  expect_equal(unname(ctrl_dims(red)["states"]), 1L)
})

test_that("arithmetic builds the expected models", {
  s <- ctrl_tf("s")
  G <- 1 / (s + 1)
  expect_equal(G$den[[1]], c(1, 1))
  H <- G^2 + 2 * G
  expect_equal(drop(ctrl_dcgain(H)), 3, tolerance = 1e-12)
  A <- ctrl_tf2ss(G)
  expect_equal(drop(ctrl_dcgain(A * A)), 1, tolerance = 1e-12)
  expect_equal(drop(ctrl_dcgain(A - A)), 0, tolerance = 1e-12)
  expect_equal(drop(ctrl_dcgain(-A)), -1, tolerance = 1e-12)
})
