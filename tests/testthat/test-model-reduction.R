test_that("balanced truncation respects the Hankel error bound", {
  G <- ctrl_tf(c(1, 11, 30), c(1, 12, 47, 72, 36))
  for (r in 1:3) {
    Gr <- ctrl_balred(G, r)
    err <- as.numeric(ctrl_norm(ctrl_tf2ss(G) - Gr, "inf"))
    expect_lte(err, attr(Gr, "error_bound") * (1 + 1e-8) + 1e-10)
    expect_true(ctrl_is_stable(Gr))
    Gm <- ctrl_balred(G, r, method = "matchdc")
    expect_equal(ctrl_dcgain(Gm), ctrl_dcgain(G), tolerance = 1e-8,
                 ignore_attr = TRUE)
  }
  hsv <- ctrl_hsvd(G)
  bal <- ctrl_balreal(G)
  expect_equal(ctrl_gram(bal, "c"), diag(hsv[seq_len(nrow(bal$A))]),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(ctrl_gram(bal, "o"), diag(hsv[seq_len(nrow(bal$A))]),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("unstable models keep their unstable part", {
  sys <- ctrl_ss(diag(c(-1, 2, -3, -10)), c(1, 1, 1, 1), c(1, 1, 1, 1))
  parts <- ctrl_stabsep(sys)
  w <- c(0.1, 1, 10)
  expect_equal(ctrl_freqresp(parts$stable + parts$unstable, w)$response,
               ctrl_freqresp(sys, w)$response, tolerance = 1e-12)
  red <- ctrl_balred(sys, 2)
  expect_equal(sum(Re(ctrl_poles(red)) > 0), 1L)
  expect_error(ctrl_balred(sys, 0), "unstable")
})

test_that("balred order 0 bound and matchdc on nonminimal models", {
  G <- ctrl_tf2ss(ctrl_tf(1, c(1, 3, 2)))
  r0 <- ctrl_balred(G, 0)
  expect_equal(attr(r0, "error_bound"), 2 * sum(ctrl_hsvd(G)))
  nm <- ctrl_ss(diag(c(-1, -2, -3)), c(1, 1, 0), c(1, 0, 1))
  r <- ctrl_balred(nm, 2, "matchdc")
  expect_equal(ctrl_dcgain(r), ctrl_dcgain(nm), tolerance = 1e-10)
})

test_that("reduction handles models with zero input-output dynamics", {
  z <- ctrl_ss(-1, 0, 1)
  expect_equal(ctrl_hsvd(z), 0)
  expect_equal(nrow(ctrl_balreal(z)$A), 0L)
  expect_equal(nrow(ctrl_balred(z, 0)$A), 0L)
})
