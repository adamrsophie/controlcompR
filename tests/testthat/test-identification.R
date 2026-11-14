arx_data <- function(Tn = 600, sd = 0.05, seed = 1) {
  set.seed(seed)
  u <- rnorm(Tn)
  e <- rnorm(Tn, sd = sd)
  y <- numeric(Tn)
  for (k in 3:Tn) y[k] <- 1.2 * y[k - 1] - 0.5 * y[k - 2] + 0.7 * u[k - 1] + 0.3 * u[k - 2] + e[k]
  ctrl_iddata(y, u)
}

test_that("ARX recovers true coefficients with valid uncertainty", {
  d <- arx_data()
  m <- ctrl_arx(d, na = 2, nb = 2, nk = 1)
  expect_equal(m$a, c(-1.2, 0.5), tolerance = 0.02)
  expect_equal(drop(m$b), c(0.7, 0.3), tolerance = 0.02)
  se <- sqrt(diag(m$vcov))
  z <- (c(m$a, m$b) - c(-1.2, 0.5, 0.7, 0.3)) / se
  expect_true(all(abs(z) < 4))
  expect_equal(m$noise_var, 0.05^2, tolerance = 0.2)
  r <- ctrl_resid(m, d)
  expect_true(r$ljung_box > 0.001)
  # transfer function equals z^-1 (0.7 + 0.3 z^-1) / (1 - 1.2 z^-1 + 0.5 z^-2)
  w <- 0.7
  zz <- exp(1i * w)
  want <- (m$b[1] * zz^-1 + m$b[2] * zz^-2) / (1 + m$a[1] * zz^-1 + m$a[2] * zz^-2)
  expect_equal(ctrl_freqresp(m$sys, w)$response[1, 1, 1], want, tolerance = 1e-10)
})

test_that("order selection and validation favor the true structure", {
  parts <- ctrl_idsplit(arx_data(1000))
  tab <- ctrl_arx_order(parts$estimation, 1:3, 1:3, validation = parts$validation)
  expect_equal(tab$na[1], 2)
  expect_equal(tab$nb[1], 2)
  cmp <- ctrl_compare(parts$validation,
                      under = ctrl_arx(parts$estimation, 1, 1, 1),
                      true = ctrl_arx(parts$estimation, 2, 2, 1))
  expect_gt(cmp$simulation_fit[2], cmp$simulation_fit[1])
  expect_gt(cmp$prediction_fit[2], 95)
})

test_that("subspace identification recovers poles and frequency response", {
  set.seed(2)
  true <- ctrl_ss(matrix(c(0.7, 0.2, -0.2, 0.7), 2), cbind(c(1, 0), c(0.3, 1)),
                  rbind(c(1, 1), c(0, 1)), dt = 1)
  Tn <- 800
  U <- matrix(rnorm(2 * Tn), Tn)
  Y <- ctrl_simulate(true, input = U, times = 0:(Tn - 1))$y[, , 1] +
    matrix(rnorm(2 * Tn, sd = 0.01), Tn)
  m <- ctrl_n4sid(ctrl_iddata(Y, U))
  expect_equal(m$order, 2L)
  expect_equal(sort_complex_test(ctrl_poles(m$sys)),
               sort_complex_test(ctrl_poles(true)), tolerance = 1e-2)
  w <- c(0.1, 1, 2)
  expect_equal(ctrl_freqresp(m$sys, w)$response, ctrl_freqresp(true, w)$response,
               tolerance = 0.02, ignore_attr = TRUE)
  fit <- ctrl_compare(ctrl_iddata(Y, U), m)
  expect_true(all(fit$simulation_fit > 95))
})

test_that("sysid models convert with the documented convention", {
  # sysid starts its regression at a different sample, so the estimates
  # agree closely but not exactly
  skip_if_not_installed("sysid")
  d <- arx_data()
  fit <- sysid::arx(sysid::idframe(output = d$y[, 1], input = d$u[, 1], Ts = 1),
                    c(2, 2, 1))
  conv <- ctrl_from_sysid(fit)
  own <- ctrl_arx(d, 2, 2, 1)
  w <- c(0.2, 1.3)
  expect_equal(ctrl_freqresp(conv$plant, w)$response,
               ctrl_freqresp(own$sys, w)$response, tolerance = 1e-2,
               ignore_attr = TRUE)
})

test_that("identification inputs are validated", {
  expect_error(ctrl_iddata(1:5, 1:4), "same number")
  d <- ctrl_iddata(rnorm(50), rep(1, 50))
  expect_error(ctrl_arx(d, 2, 3, 1), "persistently exciting")
})
