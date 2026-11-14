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

test_that("delay margin is the minimum over crossovers and hidden modes count", {
  L <- ctrl_tf(0.866 * c(1 / 1.732, 1), c(1, 0, 0)) *
    ctrl_tf(c(1, 24, 400), c(1, 0.4, 400))
  m <- ctrl_margin(L)
  expect_true(m$closed_loop_stable)
  expect_equal(m$delay_margin, min(m$all_pm$delay_margin))
  expect_lt(m$delay_margin, 0.1)
  Lh <- ctrl_ss(diag(c(-1, 1)), c(1, 0), c(1, 1))
  expect_false(ctrl_margin(Lh)$closed_loop_stable)
  expect_error(ctrl_disk_margin(Lh), "unstable")
})

test_that("gang of four uses S P for the load disturbance map", {
  P <- ctrl_ss(matrix(c(-1, 0.5, 0, -2), 2), diag(2), matrix(c(1, 0.3, 0.2, 1), 2))
  C <- ctrl_gain(matrix(c(2, 0.5, -0.4, 1), 2))
  g4 <- ctrl_gangof4(P, C)
  w <- 1.3
  Pj <- ctrl_freqresp(P, w)$response[, , 1]
  Cj <- matrix(c(2, 0.5, -0.4, 1), 2)
  PS <- ctrl_freqresp(g4$PS, w)$response[, , 1]
  expect_equal(PS, solve(diag(2) + Pj %*% Cj) %*% Pj, tolerance = 1e-10)
  P3 <- ctrl_ss(-1, matrix(1, 1, 3), 1)
  C3 <- ctrl_gain(matrix(1, 3, 1))
  expect_silent(ctrl_gangof4(P3, C3))
})

test_that("discrete margins list the Nyquist frequency once", {
  m <- ctrl_margin(ctrl_tf(0.3, c(1, 0.5), dt = 1))
  expect_equal(sum(abs(m$all_gm$frequency - pi) < 1e-9), 1)
})

test_that("H-infinity norm is not fooled by zeros at candidate frequencies", {
  A <- -diag(4)
  A[cbind(1:3, 2:4)] <- 1
  H <- ctrl_ss(A, c(0, 0, 0, 1), c(-2, 4, -3, 1))
  expect_equal(as.numeric(ctrl_norm(H, "inf")), 0.25, tolerance = 1e-6)
  expect_false(ctrl_robust_stability(ctrl_gain(0), ctrl_gain(1), 8 * H,
                                     type = "additive")$robustly_stable)
  expect_equal(as.numeric(ctrl_norm(ctrl_ss(-1, 0, 1), "inf")), 0)
})

test_that("disk phase margin follows the disk geometry", {
  dm <- ctrl_disk_margin(ctrl_tf(0.5, 1))
  expect_equal(dm$alpha, 6, tolerance = 1e-8)
  expect_equal(dm$phase_margin, 2 * atan(3) * 180 / pi, tolerance = 1e-8)
  L <- ctrl_tf(25, c(1, 10, 10, 10))
  dm <- ctrl_disk_margin(L)
  expect_equal(dm$phase_margin, 2 * atan(dm$alpha / 2) * 180 / pi, tolerance = 1e-8)
  for (sk in c(-0.5, 0.5)) {
    g <- controlcompR:::disk_gain_phase(0.8, sk)
    a <- (1 - sk) / 2; b <- (1 + sk) / 2
    d <- 0.8 * exp(1i * seq(0, 2 * pi, length.out = 20001))
    z <- (1 + a * d) / (1 - b * d)
    # largest phase on the unit circle reached by the disk boundary
    on <- abs(Mod(z) - 1) < 1e-3
    expect_equal(g$pm, max(abs(Arg(z[on]))) * 180 / pi, tolerance = 0.05)
  }
})
