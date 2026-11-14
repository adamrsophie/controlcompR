test_that("series, parallel and feedback match transfer-function algebra", {
  G1 <- ctrl_tf(1, c(1, 1))
  G2 <- ctrl_tf(c(1, 2), c(1, 3))
  S1 <- ctrl_tf2ss(G1)
  S2 <- ctrl_tf2ss(G2)
  w <- c(0.1, 1, 10)
  fr <- function(x) ctrl_freqresp(x, w)$response
  expect_equal(fr(ctrl_series(S1, S2)), fr(G1) * fr(G2), tolerance = 1e-12)
  expect_equal(fr(ctrl_parallel(S1, S2)), fr(G1) + fr(G2), tolerance = 1e-12)
  expect_equal(fr(ctrl_feedback(S1, S2)), fr(G1) / (1 + fr(G1) * fr(G2)),
               tolerance = 1e-12)
  expect_equal(fr(ctrl_feedback(G1, G2, sign = 1)), fr(G1) / (1 - fr(G1) * fr(G2)),
               tolerance = 1e-12)
  expect_equal(fr(ctrl_feedback(S1, S2)), fr(ctrl_feedback(G1, G2)), tolerance = 1e-12)
})

test_that("feedback with direct feedthrough resolves the algebraic loop", {
  G <- ctrl_ss(-1, 1, 1, D = 2)
  H <- ctrl_gain(0.5)
  cl <- ctrl_feedback(G, H)
  w <- c(0.2, 3)
  g <- ctrl_freqresp(G, w)$response
  expect_equal(ctrl_freqresp(cl, w)$response, g / (1 + 0.5 * g), tolerance = 1e-12)
  expect_error(ctrl_feedback(ctrl_gain(1), ctrl_gain(1), sign = 1), "not well posed")
})

test_that("named connections reproduce unity feedback", {
  G <- ctrl_tf(1, c(1, 1), input_names = "u", output_names = "y")
  K <- ctrl_tf(5, 1, input_names = "e", output_names = "u")
  S <- ctrl_sumblock("e = r - y")
  T <- ctrl_connect(G, K, S, inputs = "r", outputs = c("y", "u"))
  expect_equal(drop(ctrl_dcgain(T)), c(5 / 6, 5 / 6), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_error(ctrl_connect(G, K, inputs = "r", outputs = "y"), "neither")
  expect_error(ctrl_series(ctrl_ss(-1, 1, 1), ctrl_ss(-1, 1, 1, dt = 1)), "time domains")
})

test_that("numeric gains adopt the time domain of discrete models", {
  L <- ctrl_tf(0.1, c(1, -0.5), dt = 0.1)
  T <- ctrl_feedback(L, 1)
  expect_equal(T$dt, 0.1)
  expect_equal(ctrl_series(2, L)$dt, 0.1)
  expect_equal((L + 1)$dt, 0.1)
  expect_true(ctrl_margin(L)$closed_loop_stable)
  Ls <- ctrl_tf2ss(L)
  expect_equal(ctrl_feedback(Ls, 1)$dt, 0.1)
  expect_error(ctrl_series(ctrl_tf(1, c(1, 1)), L), "time domains")
})

test_that("append lets static gains adopt the time domain in any order", {
  Gd <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.1)
  expect_equal(ctrl_append(ctrl_gain(1), Gd)$dt, 0.1)
  expect_equal(ctrl_append(Gd, ctrl_gain(1))$dt, 0.1)
})

test_that("transfer function feedback checks the algebraic loop", {
  expect_error(ctrl_feedback(ctrl_tf(c(1, 1), c(1, 0)), 1, sign = 1), "algebraic loop")
  expect_error(ctrl_feedback(ctrl_tf2ss(ctrl_tf(c(1, 1), c(1, 0))), 1, sign = 1),
               "algebraic loop")
})
