# Scale, timing and validation contracts across modules.

test_that("LQ weights in tiny units give the same gain (NC05)", {
  ans <- ctrl_care(-1, 1, 1e-15, 1e-15)
  expect_equal(as.numeric(ans$K), sqrt(2) - 1, tolerance = 1e-9)
  expect_equal(as.numeric(ans$P) / 1e-15, sqrt(2) - 1, tolerance = 1e-9)
  expect_error(ctrl_care(-1, 1, 1, 0), "positive definite")
})

test_that("MPC accepts uniformly rescaled weights (MN10)", {
  p <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0)), 0.1)
  u <- vapply(c(1, 1e-13, 1e8), function(s) {
    ctrl_mpc_step(ctrl_mpc(p, 10, Q = s * diag(2), R = s, u_min = -1, u_max = 1),
                  c(1, 0))$u
  }, 0)
  expect_equal(u, rep(u[1], 3), tolerance = 1e-5)
})

test_that("MPC simulation applies the input disturbance through D (MN01)", {
  plant <- ctrl_ss(0, 0, 0, D = 1, dt = 1)
  sim <- ctrl_simulate_mpc(ctrl_mpc(plant, 2, Q = 0, Qy = 1, R = 1), 0, 3,
                           disturbance = 1)
  expect_equal(as.numeric(sim$y), as.numeric(sim$u) + 1, tolerance = 1e-8)
  # a time-varying disturbance on a MIMO plant with feedthrough
  pl <- ctrl_ss(diag(c(0.5, 0.2)), diag(2), diag(2), D = matrix(c(1, 0, 0.5, 1), 2), dt = 1)
  d <- cbind(1:4, -(1:4))
  sim <- ctrl_simulate_mpc(ctrl_mpc(pl, 3, Q = diag(2), R = diag(2)), c(1, -1), 4,
                           disturbance = d)
  for (k in 1:4) {
    expect_equal(sim$y[k, ], drop(pl$C %*% sim$x[k, ] + pl$D %*% (sim$u[k, ] + d[k, ])),
                 tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(sim$x[k + 1, ], drop(pl$A %*% sim$x[k, ] + pl$B %*% (sim$u[k, ] + d[k, ])),
                 tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("MHE estimates do not depend on a common covariance scale (MN02)", {
  plant <- ctrl_ss(1, 0, 1, dt = 1)
  est <- vapply(c(1, 1e-13, 1e8), function(s) {
    ctrl_mhe_run(ctrl_mhe(plant, 1, Qn = s, Rn = s, x0 = 0, P0 = s), 1)$estimates[1, 1]
  }, 0)
  expect_equal(est, rep(0.5, 3), tolerance = 1e-6)
})

test_that("growing MHE windows keep their dual warm start (MN09)", {
  expect_equal(shift_blocks(1:4, 2, 2, 3, 0), c(1, 2, 3, 4, 0, 0))
  expect_equal(shift_blocks(1:4, 2, 2, 2, 1), c(3, 4, 0, 0))
  expect_equal(drop_head(1:3, 0), 1:3)
})

test_that("singular system pencils are rejected (NC03)", {
  expect_error(ctrl_zeros(ctrl_ss(-1, 1, 0, 0)), "singular")
})

test_that("invariant zeros do not depend on state units (NC04)", {
  expect_equal(ctrl_zeros(ctrl_ss(-1, 1 / 2e6, 2e6, 1)), -2 + 0i, tolerance = 1e-8)
  expect_equal(ctrl_zeros(ctrl_ss(-1, 1, 1, 1)), -2 + 0i, tolerance = 1e-12)
})

test_that("H-infinity norms are invariant to state scaling (NC06)", {
  for (zeta in c(1e-4, 0.1)) {
    for (k in c(1e-6, 1e6)) {
      A0 <- matrix(c(-zeta, 1, -(1 - zeta^2), -zeta), 2, byrow = TRUE)
      Tm <- diag(c(k, 1))
      G <- ctrl_ss(solve(Tm, A0 %*% Tm), solve(Tm, c(0, 1)), matrix(c(1, 0), 1) %*% Tm)
      expect_equal(as.numeric(ctrl_norm(G, "inf")), 1 / (2 * zeta * sqrt(1 - zeta^2)),
                   tolerance = 1e-7)
    }
  }
})

test_that("state-space to transfer conversion keeps weak gains (RT02)", {
  G <- ctrl_ss(-1, 1, 1e-14)
  for (minimal in c(FALSE, TRUE)) {
    tf <- ctrl_ss2tf(G, minimal = minimal)
    expect_equal(as.complex(ctrl_freqresp(tf, 1)$response) / 1e-14, 1 / (1 + 1i),
                 tolerance = 1e-9)
  }
  expect_equal(ctrl_ss2tf(ctrl_ss(-1, 1, 0))$num[[1]], 0)
})

test_that("monic normalization rejects overflow (RT11)", {
  expect_error(ctrl_tf(1e308, c(1e-308, 1)), "overflow")
})

test_that("balanced reduction bounds dominate the error in any state units (NC01)", {
  G <- ctrl_ss(-diag(2), diag(c(1e4, 1e-4)), diag(c(1e-4, 1e4)))
  expect_equal(ctrl_hsvd(G), c(0.5, 0.5), tolerance = 1e-10)
  Gr <- ctrl_balred(G, 1)
  dc_error <- max(svd(diag(2) - ctrl_dcgain(Gr), 0, 0)$d)
  expect_gte(attr(Gr, "error_bound") + 1e-10, dc_error)
  expect_equal(ctrl_dcgain(ctrl_balred(G, 2)), diag(2), tolerance = 1e-10,
               ignore_attr = TRUE)
  # an equivalent coordinate transformation keeps the Hankel singular values
  set.seed(3)
  M <- matrix(rnorm(16), 4)
  A <- -(M %*% t(M)) - diag(4)
  S <- ctrl_ss(A, matrix(rnorm(8), 4), matrix(rnorm(8), 2))
  Q <- qr.Q(qr(matrix(rnorm(16), 4)))
  S <- ctrl_ss(t(Q) %*% A %*% Q, t(Q) %*% S$B, S$C %*% Q)
  Tm <- diag(10^c(-5, -1, 2, 6))
  S2 <- ctrl_ss(solve(Tm, S$A %*% Tm), solve(Tm, S$B), S$C %*% Tm)
  expect_equal(ctrl_hsvd(S2), ctrl_hsvd(S), tolerance = 1e-6)
  for (r in 1:3) {
    Sr <- ctrl_balred(S2, r)
    err <- ctrl_norm(S - Sr, "inf")
    expect_lte(err, attr(Sr, "error_bound") * (1 + 1e-6))
  }
  # a realization too ill-conditioned for reliable Gramians is rejected
  # instead of returning an untrustworthy bound
  Tb <- diag(10^c(-5, -1, 2, 6)) %*% qr.Q(qr(matrix(rnorm(16), 4)))
  Sb <- ctrl_ss(solve(Tb, S$A %*% Tb), solve(Tb, S$B), S$C %*% Tb)
  expect_error(suppressWarnings(ctrl_balred(Sb, 3)), "certify")
})

test_that("complex pole pairs are placed with full actuation (NC02)", {
  K <- ctrl_place(diag(c(-1, -2)), diag(2), c(-3 + 1i, -3 - 1i))
  expect_false(is.complex(K))
  ev <- eigen(diag(c(-1, -2)) - K, only.values = TRUE)$values
  expect_equal(sort(Re(ev)), c(-3, -3), tolerance = 1e-8)
  expect_equal(sort(Im(ev)), c(-1, 1), tolerance = 1e-8)
  # a repeated complex pair with two inputs
  A <- diag(c(1, 2, 3, 4))
  B <- matrix(c(1, 0, 1, 1, 0, 1, 1, 2), 4)
  K <- ctrl_place(A, B, rep(c(-1 + 1i, -1 - 1i), 2))
  ev <- eigen(A - B %*% K, only.values = TRUE)$values
  expect_equal(sort(Re(ev)), rep(-1, 4), tolerance = 1e-6)
  # observer duality
  C <- t(B)
  L <- t(ctrl_place(t(A), t(C), c(-2 + 1i, -2 - 1i, -3, -4)))
  ev <- eigen(A - L %*% C, only.values = TRUE)$values
  expect_equal(sort(Re(ev)), c(-4, -3, -2, -2), tolerance = 1e-6)
})

test_that("Nyquist encirclements match closed-loop pole counts (RT01)", {
  for (k in c(0.5, 7, 9, 1e4, 1e9)) {
    L <- ctrl_tf(k, c(1, 3, 3, 1))
    rhp <- sum(Re(polyroot(rev(c(1, 3, 3, 1 + k)))) > 0)
    expect_equal(ctrl_nyquist(L)$encirclements, rhp)
  }
  L <- ctrl_c2d(ctrl_tf(20, c(1, 3, 3, 1)), 0.1)
  expect_equal(ctrl_nyquist(L)$encirclements,
               sum(Mod(ctrl_poles(ctrl_feedback(L, 1))) > 1))
  # open-loop unstable plant: Z = N + P
  L <- ctrl_tf(3, c(1, -1))
  ny <- ctrl_nyquist(L)
  expect_equal(ny$encirclements + ny$open_loop_unstable, 0)
  # 1 + D = 0 makes the contour singular at infinity
  expect_true(is.na(ctrl_nyquist(ctrl_tf(c(-1, 0), c(1, 1)))$encirclements))
})

test_that("robust performance checks the frequency endpoints (ERI01)", {
  r <- ctrl_robust_performance(ctrl_gain(0), ctrl_gain(0),
                               ctrl_tf(c(1.001, 0), c(1, 1)), 0)
  expect_false(r$robust_performance)
  expect_false(r$certified)
  expect_equal(r$peak, 1.001, tolerance = 1e-9)
  # discrete loop whose largest weighted sensitivity is at DC
  Wd <- ctrl_tf(1.01, c(1, -0.5), dt = 1)  # |Wd| peaks at z = 1 with 2.02
  r <- ctrl_robust_performance(ctrl_gain(0, dt = 1), ctrl_gain(0, dt = 1), Wd, 0)
  expect_equal(r$peak, 2.02, tolerance = 1e-9)
  expect_false(r$robust_performance)
  # a comfortable pass is certified by the rigorous upper bound
  G <- ctrl_tf(1, c(1, 1))
  K <- ctrl_gain(1)
  r <- ctrl_robust_performance(G, K, 0.2, 0.2)
  expect_true(r$robust_performance)
  expect_true(r$certified)
  expect_gte(r$upper_bound, r$peak)
})

test_that("trajectory optimization reports the soft penalty in its cost (MN03)", {
  sys <- ctrl_nlsys(function(k, x, u, params) x + u, nstates = 1, ninputs = 1, dt = 1)
  tr <- ctrl_nl_trajopt(sys, 1, horizon = 1, Q = 0, R = 1, P = 0, x_max = 0, penalty = 1)
  expect_equal(as.numeric(tr$u), -0.5, tolerance = 1e-5)
  expect_equal(tr$objective, 0.5, tolerance = 1e-6)
  expect_equal(unname(tr$cost[["penalty"]]), 0.25, tolerance = 1e-6)
})

test_that("NMPC hold fallback respects the input bounds (MN04)", {
  sys <- ctrl_nlsys(function(k, x, u, params) x + u, nstates = 1, ninputs = 1, dt = 1)
  # maxit = 0 makes every optimization report failure
  nm <- ctrl_nmpc(sys, 5, Q = 1, R = 1, u_min = 1, u_max = 2, maxit = 0)
  sm <- suppressWarnings(ctrl_simulate_nmpc(nm, -10, 3, on_failure = "hold"))
  expect_true(all(sm$u >= 1 & sm$u <= 2))
  expect_true(all(sm$log$fallback != "none"))
  # the default held input is the admissible point closest to zero
  expect_equal(as.numeric(sm$u[1, ]), 1)
  sm <- suppressWarnings(ctrl_simulate_nmpc(nm, -10, 2, on_failure = "hold", u0 = 1.5))
  expect_equal(as.numeric(sm$u), c(1.5, 1.5))
  expect_error(ctrl_simulate_nmpc(nm, -10, 2, on_failure = "hold", u0 = 0), "bounds")
})

test_that("nonlinear filters validate covariances and data (MN05)", {
  sys <- ctrl_nlsys(function(k, x, u, params) x, nstates = 1, ninputs = 0, dt = 1)
  expect_error(ctrl_ekf(sys, 1, Qn = 0, Rn = 1, x0 = 0, P0 = -0.5), "P0")
  expect_error(ctrl_ukf(sys, 1, Qn = 0, Rn = 1, x0 = 0, P0 = -0.5), "P0")
  expect_error(ctrl_ekf(sys, c(0, Inf), Qn = 0, Rn = 1), "infinite")
  expect_error(ctrl_ukf(sys, 1, Qn = 0, Rn = 1, alpha = 0), "alpha")
})

test_that("RK4 discretization needs a positive integer number of substeps (MN06)", {
  sys <- ctrl_nlsys(function(t, x, u, p) 1, nstates = 1, ninputs = 0)
  expect_error(ctrl_nl_discretize(sys, 1, "rk4", substeps = 0), "substeps")
  expect_error(ctrl_nl_discretize(sys, 1, "rk4", substeps = 2.5), "substeps")
})

test_that("the NMPC simulator rejects incompatible plants (MN07)", {
  discrete <- ctrl_nlsys(function(k, x, u, p) 0.9 * x + u, nstates = 1, ninputs = 1, dt = 1)
  continuous <- ctrl_nlsys(function(t, x, u, p) -x + u, nstates = 1, ninputs = 1)
  nm <- ctrl_nmpc(discrete, 2, Q = 1, R = 1)
  expect_error(ctrl_simulate_nmpc(nm, 1, 2, plant = continuous), "discrete")
  other <- ctrl_nlsys(function(k, x, u, p) 0.9 * x + u, nstates = 1, ninputs = 1, dt = 0.5)
  expect_error(ctrl_simulate_nmpc(nm, 1, 2, plant = other), "sample interval")
})

test_that("an unconverged equilibrium search is not classified (MN08)", {
  sys <- ctrl_nlsys(function(t, x, u, p) exp(-x), nstates = 1, ninputs = 0)
  eq <- ctrl_equilibrium(sys, 0, maxit = 1)
  expect_false(eq$converged)
  expect_match(eq$stability, "not assessed")
})

test_that("LQG simulation reports the controller's own estimate (ERI07, ERI08)", {
  design <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, -1, -0.5), 2, byrow = TRUE), c(0, 1),
                             c(1, 0)), 0.1)
  lqg <- ctrl_lqg(design, Q = diag(2), R = 1, Qn = 0.01 * diag(2), Rn = 0.1,
                  form = "current")
  actual <- design
  actual$C <- matrix(c(1.3, 0.2), 1)
  sim <- ctrl_simulate_lqg(actual, lqg, steps = 30, x0 = c(1, 0), seed = 1)
  K <- lqg$lqr$K
  u_from_xhat <- -drop(sim$xhat[1:30, , 1] %*% t(K))
  expect_equal(u_from_xhat, sim$u[, 1, 1], tolerance = 1e-8)
  slow <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, -1, -0.5), 2, byrow = TRUE), c(0, 1),
                           c(1, 0)), 1)
  expect_error(ctrl_simulate_lqg(slow, lqg, steps = 5, x0 = c(1, 0)), "sample interval")
})

test_that("infinite measurements are not treated as missing (ERI09)", {
  sys <- ctrl_ss(1, 1, 1, dt = 1)
  expect_error(ctrl_kalman_filter(sys, c(0, Inf, 0), Qn = 1, Rn = 1), "infinite")
  expect_silent(ctrl_kalman_filter(sys, c(0, NA, 0), Qn = 1, Rn = 1))
})

test_that("nonlinear filter results are not passed to the linear smoother (ERI10)", {
  sys <- ctrl_nlsys(function(k, x, u, p) x, nstates = 1, ninputs = 0, dt = 1)
  f <- ctrl_ekf(sys, c(1, 2), Qn = 1, Rn = 1)
  expect_s3_class(f, "ctrl_nl_filter")
  expect_error(ctrl_kalman_smoother(f), "only implemented for the linear")
})

test_that("mixed-sensitivity weights keep their time domain (ERI06)", {
  G <- ctrl_tf(1, c(1, 1))
  expect_error(ctrl_augw(G, W1 = ctrl_tf(1, c(1, -0.5), dt = 1)), "discrete time weight")
  P <- ctrl_augw(G, W1 = ctrl_gain(2, dt = 1))
  expect_equal(P$dt, 0)
})

test_that("ARX order selection compares models on common responses (ERI02)", {
  set.seed(194)
  y <- rnorm(200)
  u <- rnorm(200)
  a <- ctrl_arx_order(ctrl_iddata(y, u), na = 1:6, nb = 1)
  b <- ctrl_arx_order(ctrl_iddata(1000 * y, u), na = 1:6, nb = 1)
  expect_identical(a$na, b$na)
  expect_true(all(a$nobs == a$nobs[1]))
  expect_error(ctrl_arx(ctrl_iddata(y, u), 2, 1, start = 1), "start")
})

test_that("N4SID order selection and identifiability checks (ERI03, ERI04)", {
  set.seed(7)
  plant <- ctrl_ss(diag(c(0.2, 0.8)), matrix(c(1, 1)), matrix(c(1, 1), 1), dt = 1)
  u <- rnorm(300)
  y <- ctrl_simulate(plant, u, 0:299)$y[, 1, 1]
  fit <- ctrl_n4sid(ctrl_iddata(y, u), horizon = 3)
  expect_equal(fit$order, 2L)
  expect_equal(sort(Mod(ctrl_poles(fit$sys))), c(0.2, 0.8), tolerance = 1e-6)
  expect_error(ctrl_n4sid(ctrl_iddata(y, u), order = 3, horizon = 3), "horizon")
})

test_that("residual analysis uses the lags the data support (ERI05)", {
  set.seed(3)
  u <- rnorm(30)
  y <- 0.5 * c(0, u[-30]) + rnorm(30, sd = 0.1)
  d <- ctrl_iddata(y, u)
  m <- ctrl_arx(d, 1, 1)
  r <- ctrl_resid(m, d, lags = 100)
  expect_equal(r$lags, nrow(as.matrix(r$residuals)) - 1L)
  expect_equal(max(r$acf$lag), r$lags)
  expect_error(ctrl_resid(m, d, lags = 0), "lags")
})

test_that("standard responses start at time zero (RT06)", {
  G <- ctrl_tf(1, c(1, 1))
  st <- ctrl_step(G, times = 1:3)
  expect_equal(st$time, 1:3)
  expect_equal(st$y[, 1, 1], 1 - exp(-(1:3)), tolerance = 1e-8)
  expect_equal(ctrl_impulse(G, times = 1:3)$y[, 1, 1], exp(-(1:3)), tolerance = 1e-8)
  expect_equal(ctrl_initial(G, 1, times = 2:3)$y[, 1, 1],
               ctrl_initial(G, 1, times = 0:3)$y[3:4, 1, 1])
  Gd <- ctrl_c2d(G, 0.5)
  expect_equal(ctrl_step(Gd, times = c(1, 1.5))$y[, 1, 1],
               ctrl_step(Gd, times = seq(0, 1.5, 0.5))$y[3:4, 1, 1])
  expect_error(ctrl_step(G, times = c(-1, 0)), "at or after 0")
})

test_that("an exactly known final value does not imply settling (RT07)", {
  G <- ctrl_tf(1, c(1, 1))
  info <- ctrl_step_info(ctrl_step(G, times = seq(0, 0.02, by = 0.001)), sys = G)
  expect_true(info$final_value_known)
  expect_false(info$settled)
  expect_true(is.na(info$settling_time))
  expect_true(ctrl_step_info(G)$settled)
})

test_that("FOPDT fitting normalizes by the step amplitude (RT08)", {
  G <- ctrl_tf(2, c(1, 1))
  st <- ctrl_step(G, times = seq(0, 10, by = 0.01), amplitude = 5)
  expect_equal(unname(ctrl_fit_fopdt(st)["K"]), 2, tolerance = 1e-3)
})

test_that("PID initialization modes (RT09)", {
  pid <- ctrl_pid(kp = 0, ki = 1, dt = 1, integral_method = "tustin")
  expect_equal(ctrl_pid_update(pid, NULL, r = 1, y = 0)$u, 1)
  z <- ctrl_pid_update(pid, NULL, r = 1, y = 0, init = "zero")$u
  tf <- ctrl_pid_tf(pid)
  first <- ctrl_simulate(tf, rep(1, 3), 0:2)$y[1, 1, 1]
  expect_equal(z, first)
  expect_equal(z, 0.5)
})

test_that("an unbounded affine PID loop with negative gain is solved (RT10)", {
  plant <- ctrl_gain(-2, dt = 1)
  pid <- ctrl_pid(kp = 1, dt = 1)
  sim <- ctrl_simulate_pid(plant, pid, reference = 1, steps = 2)
  expect_equal(sim$u[1], -1)
  expect_equal(sim$y[1], 2)
  pidb <- ctrl_pid(kp = 1, dt = 1, u_min = -5, u_max = 5)
  expect_error(ctrl_simulate_pid(plant, pidb, reference = 1, steps = 2), "positive")
})

test_that("input callbacks must return a full sample every time (RT12)", {
  G <- ctrl_ss(-diag(2), diag(2), diag(2))
  f <- function(t) if (t == 0) c(1, 2) else 3
  expect_error(ctrl_simulate(G, f, seq(0, 1, by = 0.5)), "returned 1")
})

test_that("performance metrics need increasing time (RT13) and keep the sampled convention (RT05)", {
  expect_error(ctrl_performance(c(2, 1, 0), c(0, 0, 0), 1), "increasing")
  sims <- list(a = function() list(time = 0:2, y = c(0, 0, 0)))
  out <- ctrl_compare_controllers(sims, reference = 1, discrete = TRUE)
  expect_equal(out$IAE, 3)
  expect_equal(out$convention, "sampled")
  sims$b <- function() list(time = 0:2, y = c(0, 0, 0), discrete = TRUE)
  out <- ctrl_compare_controllers(sims, reference = 1)
  expect_equal(out$IAE, c(2, 3))
})

test_that("disk margins follow the geometry for negative skew (NC07)", {
  dm <- ctrl_disk_margin(ctrl_tf(-0.25, 1), skew = -2)
  expect_equal(as.numeric(c(dm$gain_margin, dm$phase_margin)), c(0, 2.5, 180),
               tolerance = 1e-8)
  # skew -1: the disk 1 + delta, |delta| <= alpha
  d <- ctrl_disk_margin(ctrl_tf(4, c(1, 2, 1)), skew = -1)
  expect_equal(d$gain_margin, c(1 - d$alpha, 1 + d$alpha), tolerance = 1e-6)
})

test_that("improper feedback factors are handled as rational algebra (RT03)", {
  H <- ctrl_feedback(-ctrl_tf("s"), 1)
  expect_equal(H$num[[1]], c(1, 0) * sign(H$den[[1]][1]), tolerance = 1e-12)
})

test_that("numeric scalars scale MIMO models (RT04)", {
  G <- ctrl_ss(-diag(2), diag(2), matrix(1:6, 3, 2))
  expect_equal(ctrl_dcgain(2 * G), 2 * ctrl_dcgain(G), ignore_attr = TRUE)
  expect_equal(ctrl_dcgain(G * 2), 2 * ctrl_dcgain(G), ignore_attr = TRUE)
  expect_equal(ctrl_dcgain(G / 2), ctrl_dcgain(G) / 2, ignore_attr = TRUE)
  expect_equal((2 * G)$output_names, G$output_names)
  expect_error(G / 0, "zero")
})
