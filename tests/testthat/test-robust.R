lqg_plant <- function() {
  A <- matrix(c(0, 1, 0, 0), 2, byrow = TRUE)
  ctrl_ss(A, cbind(diag(2), 0, c(0, 1)),
          rbind(c(1, 0), c(0, 0), c(1, 0)),
          rbind(c(0, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0)))
}

test_that("H2 synthesis reproduces continuous LQG", {
  P <- lqg_plant()
  h2 <- ctrl_h2syn(P, 1, 1)
  A <- P$A
  plant <- ctrl_ss(A, c(0, 1), c(1, 0))
  lqg <- ctrl_lqg(plant, Q = diag(c(1, 0)), R = 1, Qn = diag(2), Rn = 1)
  # LQG controller maps y to u with u = K y in the LFT convention
  fr1 <- ctrl_freqresp(h2$K, c(0.1, 1, 10))$response
  fr2 <- ctrl_freqresp(lqg$controller, c(0.1, 1, 10))$response
  expect_equal(fr1, fr2, tolerance = 1e-7, ignore_attr = TRUE)
  expect_true(ctrl_is_stable(h2$closed_loop))
  # local optimality: perturbing the controller does not reduce the norm
  set.seed(1)
  k <- as_ss(h2$K)
  for (i in 1:5) {
    kp <- ctrl_ss(k$A, k$B + 1e-3 * rnorm(length(k$B)), k$C, k$D)
    clp <- ctrl_lft(P, kp, 1, 1)
    expect_gte(ctrl_norm(clp, "2"), h2$norm - 1e-9)
  }
})

test_that("H-infinity synthesis is certified and brackets the optimum", {
  P <- lqg_plant()
  hinf <- ctrl_hinfsyn(P, 1, 1, tol = 1e-4)
  expect_lt(hinf$norm, hinf$gamma)
  expect_true(ctrl_is_stable(hinf$closed_loop))
  expect_error(ctrl_hinfsyn(P, 1, 1, gamma = 0.95 * hinf$gamma), "not achievable")
  h2 <- ctrl_h2syn(P, 1, 1)
  expect_lte(hinf$norm, ctrl_norm(h2$closed_loop, "inf") + 1e-8)
  # large gamma: central controller tends to the H2 controller
  big <- ctrl_hinfsyn(P, 1, 1, gamma = 1e4)
  fr1 <- ctrl_freqresp(big$K, c(0.3, 3))$response
  fr2 <- ctrl_freqresp(h2$K, c(0.3, 3))$response
  expect_equal(fr1, fr2, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("H-infinity handles scaling, cross terms and D22", {
  set.seed(8)
  A <- matrix(rnorm(9), 3)
  B1 <- matrix(rnorm(6), 3)
  B2 <- matrix(rnorm(3), 3)
  C1 <- matrix(rnorm(6), 2)
  C2 <- matrix(rnorm(3), 1)
  D12 <- matrix(c(0.3, 2), 2)
  D21 <- matrix(c(0.5, 1.5), 1)
  D22 <- matrix(0.4)
  P <- ctrl_ss(A, cbind(B1, B2), rbind(C1, C2),
               rbind(cbind(matrix(0, 2, 2), D12), cbind(D21, D22)))
  h2 <- ctrl_h2syn(P, 1, 1)
  expect_true(ctrl_is_stable(h2$closed_loop))
  hinf <- ctrl_hinfsyn(P, 1, 1)
  expect_lt(hinf$norm, hinf$gamma)
  expect_lte(hinf$norm, ctrl_norm(h2$closed_loop, "inf") + 1e-8)
})

test_that("unsupported synthesis problems are rejected", {
  P <- lqg_plant()
  Pd <- ctrl_c2d(P, 0.1)
  expect_error(ctrl_hinfsyn(Pd, 1, 1), "continuous")
  P2 <- P
  P2$D[1, 1] <- 1
  expect_error(ctrl_hinfsyn(P2, 1, 1), "D11")
  P3 <- P
  P3$D[2, 4] <- 0
  expect_error(ctrl_h2syn(P3, 1, 1), "D12")
})

test_that("mixed sensitivity design is certified and shapes S", {
  G <- ctrl_tf(1, c(1, 3, 2))
  W1 <- ctrl_tf(1, c(1, 0.01))
  P <- ctrl_augw(G, W1 = W1, W2 = 0.1)
  syn <- ctrl_hinfsyn(P, 1, 1)
  expect_lt(syn$norm, syn$gamma)
  # closed loop with u = K e equals negative feedback with -K... check S bound
  K <- syn$K
  L <- ctrl_series(K, G)
  S <- ctrl_feedback(ctrl_gain(1), L)
  w <- 10^seq(-3, 1, length.out = 50)
  sw <- Mod(ctrl_freqresp(S, w)$response[1, 1, ] *
              ctrl_freqresp(W1, w)$response[1, 1, ])
  expect_lt(max(sw), syn$gamma + 1e-6)
})

test_that("small-gain robust stability test", {
  G <- ctrl_tf(1, c(1, 1))
  K <- ctrl_tf(10, 1)
  rs <- ctrl_robust_stability(G, K, 0.5)
  # T = 10 / (s + 11): peak 10/11
  expect_equal(rs$peak, 0.5 * 10 / 11, tolerance = 1e-6)
  expect_true(rs$robustly_stable)
  rs2 <- ctrl_robust_stability(G, K, 2)
  expect_false(rs2$robustly_stable)
  ra <- ctrl_robust_stability(G, K, 0.05, type = "additive")
  # KS = 10 (s + 1)/(s + 11): peak 10 at high frequency
  expect_equal(ra$peak, 0.05 * 10, tolerance = 1e-6)
  unstable <- ctrl_robust_stability(G, ctrl_tf(-5, 1), 0.1)
  expect_false(unstable$nominally_stable)
})

test_that("hinfsyn validates the tolerance", {
  P <- ctrl_ss(-1, cbind(1, 1), rbind(1, 1), rbind(c(0, 1), c(1, 0)))
  expect_error(ctrl_hinfsyn(P, 1, 1, tol = -1), "tol")
})

test_that("robust performance matches the SISO condition", {
  G <- ctrl_tf(1, c(1, 1))
  K <- ctrl_tf(c(2, 1), c(1, 0))
  Wp <- ctrl_tf(c(0.5, 0.5), c(1, 0.01))
  Wu <- ctrl_tf(c(1, 0.2), c(0.1, 1))
  rp <- ctrl_robust_performance(G, K, Wp, Wu)
  w <- 10^seq(-4, 4, length.out = 200001)
  s <- 1i * w
  Gw <- 1 / (s + 1); Kw <- (2 * s + 1) / s
  S <- 1 / (1 + Gw * Kw); T <- 1 - S
  ref <- max(Mod((0.5 * s + 0.5) / (s + 0.01)) * Mod(S) + Mod((s + 0.2) / (0.1 * s + 1)) * Mod(T))
  expect_equal(rp$peak, ref, tolerance = 1e-6)
  expect_true(rp$criterion_exact)
  expect_equal(rp$robust_performance, ref < 1)
  expect_error(ctrl_robust_performance(G, K, ctrl_tf2ss(ctrl_append(Wp, Wp)), Wu), "scalar")
})
