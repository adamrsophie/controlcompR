test_that("classical margins match analytic values", {
  L <- ctrl_tf(10, c(1, 3, 3, 1))  # 10 / (s + 1)^3
  m <- ctrl_margin(L)
  expect_equal(m$phase_crossover, sqrt(3), tolerance = 1e-10)
  expect_equal(m$gain_margin, 0.8, tolerance = 1e-10)
  wgc <- sqrt(10^(2 / 3) - 1)
  expect_equal(m$gain_crossover, wgc, tolerance = 1e-10)
  expect_equal(m$phase_margin, 180 - 3 * atan(wgc) * 180 / pi, tolerance = 1e-8)
  expect_false(m$closed_loop_stable)
  m2 <- ctrl_margin(ctrl_tf(2, c(1, 3, 3, 1)))
  expect_equal(m2$gain_margin, 4, tolerance = 1e-10)
  expect_true(m2$closed_loop_stable)
  expect_equal(ctrl_margin(ctrl_tf2ss(ctrl_tf(2, c(1, 3, 3, 1))))$gain_margin, 4,
               tolerance = 1e-8)
})

test_that("discrete margins agree with a dense frequency grid", {
  L <- ctrl_c2d(ctrl_tf(2, c(1, 3, 3, 1)), 0.1)
  m <- ctrl_margin(L)
  w <- seq(1e-3, pi / 0.1, length.out = 200001)
  g <- ctrl_freqresp(L, w)$response[1, 1, ]
  k <- which.min(abs(Im(g)) + (Re(g) > 0) * 1e6)
  expect_equal(m$gain_margin, 1 / Mod(g[k]), tolerance = 1e-4)
  kg <- which.min(abs(Mod(g) - 1))
  expect_equal(m$gain_crossover, w[kg], tolerance = 1e-3)
})

test_that("H2 and H-infinity norms", {
  expect_equal(ctrl_norm(ctrl_tf(1, c(1, 1)), "2"), sqrt(0.5), tolerance = 1e-12)
  z <- 0.1
  G <- ctrl_tf(1, c(1, 2 * z, 1))
  expect_equal(as.numeric(ctrl_norm(G, "inf")), 1 / (2 * z * sqrt(1 - z^2)),
               tolerance = 1e-8)
  expect_equal(attr(ctrl_norm(G, "inf"), "frequency"), sqrt(1 - 2 * z^2),
               tolerance = 1e-4)
  sys <- random_stable(5, 2, 3, seed = 11)
  w <- 10^seq(-3, 3, length.out = 20000)
  fr <- ctrl_freqresp(sys, w)$response
  grid_max <- max(apply(fr, 3, function(M) max(svd(M)$d)))
  hn <- as.numeric(ctrl_norm(sys, "inf"))
  expect_gte(hn, grid_max * (1 - 1e-9))
  expect_lt(hn, grid_max * 1.01)
  d <- ctrl_c2d(G, 0.2)
  wd <- seq(1e-4, pi / 0.2, length.out = 50000)
  gd <- max(Mod(ctrl_freqresp(d, wd)$response))
  expect_equal(as.numeric(ctrl_norm(d, "inf")), gd, tolerance = 1e-4)
  # discrete H2 norm equals the root sum of squared impulse response samples
  dd <- ctrl_ss(0.5, 1, 1, dt = 1)
  expect_equal(ctrl_norm(dd, "2"), sqrt(1 / (1 - 0.25)), tolerance = 1e-12)
  expect_equal(ctrl_norm(ctrl_tf(1, c(1, -1)), "inf"), Inf)
})

test_that("disk margins match the sensitivity peak on a grid", {
  L <- ctrl_tf(25, c(1, 10, 10, 10))
  dm <- ctrl_disk_margin(L)
  w <- 10^seq(-3, 3, length.out = 50000)
  g <- ctrl_freqresp(L, w)$response[1, 1, ]
  S <- 1 / (1 + g)
  expect_equal(dm$alpha, 1 / max(Mod(S - 0.5)), tolerance = 1e-5)
  expect_equal(dm$phase_margin, 2 * atan(dm$alpha / 2) * 180 / pi, tolerance = 1e-10)
  expect_lt(dm$gain_margin[1], 1)
  expect_gt(dm$gain_margin[2], 1)
})

test_that("Bode, Nyquist, bandwidth and root locus", {
  G <- ctrl_tf(1, c(1, 1))
  b <- ctrl_bode(G, c(1))
  expect_equal(b$data$magnitude_db, -10 * log10(2), tolerance = 1e-12)
  expect_equal(b$data$phase, -45, tolerance = 1e-12)
  expect_equal(ctrl_bandwidth(G), sqrt(10^0.3 - 1), tolerance = 1e-8)
  ny <- ctrl_nyquist(ctrl_tf(10, c(1, 3, 3, 1)))
  expect_equal(ny$encirclements, 2L)
  ny2 <- ctrl_nyquist(ctrl_tf(2, c(1, 3, 3, 1)))
  expect_equal(ny2$encirclements, 0L)
  rl <- ctrl_rlocus(ctrl_tf(1, c(1, 3, 2)), gains = c(0, 0.25, 10))
  expect_equal(sort(Re(rl$roots[2, ])), c(-1.5, -1.5), tolerance = 1e-6)
  pz <- ctrl_pzmap(G)
  expect_s3_class(plot(pz), "ggplot")
  expect_s3_class(plot(b), "ggplot")
  expect_s3_class(plot(ny), "ggplot")
  expect_s3_class(plot(rl), "ggplot")
  expect_s3_class(plot(ctrl_nichols(G)), "ggplot")
  expect_s3_class(plot(ctrl_step(G)), "ggplot")
})
