# Scale, rotation and range contracts: each test asserts an analytic result
# or an explicit rejection of invalid input.


test_that("NC01 rotated equivalent realization retains its Hankel values", {
  U <- matrix(c(1,1,-1,1),2)/sqrt(2)
  B <- U%*%diag(c(1e4,1e-4))
  C <- diag(c(1e-4,1e4))%*%t(U)
  G <- ctrl_ss(-diag(2),B,C)
  expect_equal(as.numeric(ctrl_hsvd(G)),svd(C%*%B)$d/2,tolerance=1e-6)
})

test_that("NC01 balanced realization preserves nonzero transfer after rotation", {
  U <- matrix(c(1,1,-1,1),2)/sqrt(2)
  B <- U%*%diag(c(1e4,1e-4))
  C <- diag(c(1e-4,1e4))%*%t(U)
  G <- ctrl_ss(-diag(2),B,C)
  expect_equal(unname(ctrl_dcgain(ctrl_balreal(G))),unname(C%*%B),tolerance=1e-6)
})

test_that("A2-NC01 stability is invariant to an exact diagonal state rescaling", {
  A <- matrix(c(-1,1e13,0,-2),2,byrow=TRUE)
  G <- ctrl_ss(A,c(0,1e-13),c(1,0))
  expect_true(ctrl_is_stable(G))
})

test_that("A2-NC01 a uniformly slow stable scalar still has its finite norm", {
  expect_equal(as.numeric(ctrl_norm(ctrl_ss(-1e-14,1e-14,1),"inf")),1,tolerance=1e-8)
})

test_that("A2-NC02 full normal rank is not disproved by two zero locations", {
  zeros <- c(.37+1.1i,.37-1.1i,-1.3+.71i,-1.3-.71i)
  cp <- 1+0i
  for(r in zeros+1)cp <- c(cp,0)-r*c(0,cp)
  A <- -diag(4); A[cbind(1:3,2:4)] <- 1
  B <- matrix(0,4,2); B[4,2] <- 1
  C <- rbind(rep(0,4),.01*rev(Re(cp[-1])))
  G <- ctrl_ss(A,B,C,diag(c(1,.01)))
  got <- ctrl_zeros(G)
  expect_equal(sort(Re(got)),sort(Re(zeros)),tolerance=1e-7)
  expect_equal(sort(Im(got)),sort(Im(zeros)),tolerance=1e-7)
})

test_that("A2-NC03 scalar CARE is invariant to input-unit conversion", {
  k <- 1e8
  expect_equal(as.numeric(ctrl_care(-1,k,1,k^2)$P),sqrt(2)-1,tolerance=1e-6)
})

test_that("A2-NC03 scalar DARE is invariant to input-unit conversion", {
  k <- 1e8
  expect_equal(as.numeric(ctrl_dare(1,k,1,k^2)$P),(1+sqrt(5))/2,tolerance=1e-6)
})

test_that("A2-ERI01 RTS posterior is independent of unused-state units", {
  G <- ctrl_ss(diag(2), matrix(0, 2, 1), matrix(c(0, 1), 1), dt = 1)
  f <- ctrl_kalman_filter(G, c(NA_real_, 1), Qn = matrix(0, 2, 2),
                          Rn = 0.01, P0 = diag(c(1e14, 1)))
  sm <- ctrl_kalman_smoother(f)
  # x2 is constant. One N(x2, .01) observation with prior N(0, 1)
  # gives mean 1/1.01 and variance .01/1.01 at both sample times.
  expect_equal(as.numeric(sm$smoothed[, 2]), rep(1 / 1.01, 2),
                         tolerance = 1e-10)
  expect_equal(as.numeric(sm$P_smoothed[2, 2, ]), rep(0.01 / 1.01, 2),
                         tolerance = 1e-10)
})

test_that("A2-ERI02 static loops adopt a discrete dynamic weight's clock", {
  W <- ctrl_tf(0.3, c(1, -0.8), dt = 1)
  r <- ctrl_robust_performance(ctrl_gain(0), ctrl_gain(0), W, 0, n = 80)
  # S = 1, T = 0. On the unit circle the peak is at z=1.
  expect_equal(r$peak, 0.3 / (1 - 0.8), tolerance = 1e-9)
  expect_false(r$robust_performance)
  expect_false(r$nominal_performance)
  # This assertion already passes: the separate certificate protects the case.
  expect_false(r$certified)
})

test_that("A2-ERI03 automatic Hinf synthesis handles an attained zero optimum", {
  # z=u; xdot=-x+u; y=x+w. K=0 is internally stable and gives Tzw=0.
  P <- ctrl_ss(-1, cbind(0, 1), rbind(0, 1), rbind(c(0, 1), c(1, 0)))
  h <- ctrl_hinfsyn(P, 1, 1)
  expect_true(ctrl_is_stable(h$closed_loop))
  expect_equal(as.numeric(h$norm), 0, tolerance = 1e-12)
})

test_that("A2-ERI04 synthesis rank preconditions use relative block scales", {
  k <- 1e11
  P <- ctrl_ss(-1, cbind(1, 0, k), rbind(1, 0, 1),
               rbind(c(0, 0, 0), c(0, 0, k), c(0, 1, 0)))
  pp <- ctrl_plant_partition(P, 1, 1)
  # D12=[0,k]' and D21=[0,1] each have rank one and condition ratio one.
  # A=-1 is stable, B2=k is controllable, and C2=1 is observable.
  expect_equal(as.numeric(svd(pp$D12, 0, 0)$d), k)
  expect_equal(as.numeric(svd(pp$D21, 0, 0)$d), 1)
  # Isolate this reported wrapper defect from the separate input-unit
  # limitation in the downstream CARE solver (owned by the core review).
  expect_no_error(controlcompR:::synthesis_checks(pp, "H2"))
})

test_that("A2-ERI05 a nonfinite Kalman prior mean is rejected", {
  G <- ctrl_ss(0.8, 1, 1, dt = 1)
  for (x0 in c(NA_real_, Inf, NaN)) {
    expect_error(ctrl_kalman_filter(G, c(0, 1, 0), Qn = 0.1,
                                             Rn = 0.2, x0 = x0))
  }
})

test_that("UKF remains exact for an anisotropic linear prior", {
  sys <- ctrl_nlsys(function(k,x,u,p) x,function(k,x,u,p) x[2],
                    nstates=2,ninputs=0,dt=1)
  y <- c(1e-5,0)
  r <- ctrl_ukf(sys,y,Qn=0,Rn=1e-12,x0=c(0,0),P0=diag(c(1,1e-10)))
  pv <- 1/(1e10+(1:2)/1e-12)
  xm <- pv*cumsum(y)/1e-12
  # Normalize before comparison; absolute tolerances would hide covariance errors.
  expect_equal(as.numeric(r$P_filtered[2,2,]/pv),c(1,1),tolerance=1e-9)
  expect_equal(as.numeric(r$filtered[,2]/1e-5),xm/1e-5,tolerance=1e-9)
})

test_that("MHE retains an accepted positive covariance direction", {
  sys <- ctrl_ss(diag(2),matrix(0,2,1),c(0,1),dt=1)
  P0 <- diag(c(1,1e-13))
  mhe <- ctrl_mhe(sys,1,Qn=P0,Rn=1e-13,x0=c(0,0),P0=P0)
  r <- ctrl_mhe_run(mhe,1e-7)
  expect_equal(unname(r$estimates[1,2]/1e-7),0.5,tolerance=1e-7)
})

test_that("NMPC rejects negative state-violation penalties", {
  sys <- ctrl_nlsys(function(k,x,u,p) x+u,nstates=1,ninputs=1,dt=1)
  expect_error(ctrl_nmpc(sys,1,Q=0,R=1,P=0,u_min=-2,u_max=2,
                        x_max=0,penalty=-10))
})

test_that("matched discretization preserves finite DC gain across pole scales", {
  g <- ctrl_tf(1,c(1,1e8+1e-6,100))
  gd <- ctrl_c2d(g,1e7,method="matched")
  expect_equal(Re(drop(ctrl_freqresp(gd,0)$response))/0.01,1,tolerance=1e-10)
})

test_that("TF realization does not underflow a nonzero transfer to zero", {
  s <- ctrl_tf2ss(ctrl_tf(1,c(1,1e100,0,0,0)),minimal=FALSE)
  expect_true(any(s$C != 0))
})

test_that("TF realization avoids overflow from optional frequency scaling", {
  # Finite, valid polynomial; old unscaled realization returns 1 + 1e-100i at j.
  expect_no_error(ctrl_tf2ss(ctrl_tf(1,c(1,1e-100,0,0,0)),minimal=FALSE))
})

test_that("scalar arithmetic preserves the finite-model invariant", {
  expect_error(2 * ctrl_tf(1e308,1))
  expect_error(2 * ctrl_ss(-1,1,1e308,0))
})

test_that("PID simulation obeys its own plant feedthrough equation", {
  sim <- ctrl_simulate_pid(ctrl_gain(.5,dt=1),ctrl_pid(kp=1,dt=1),
                           reference=1e16,steps=2)
  expect_equal(sim$u/1e16,rep(2/3,2),tolerance=1e-12)
  expect_equal(sim$y_true/1e16,.5*sim$u/1e16,tolerance=1e-12)
})

test_that("Riccati solutions are invariant to input, cost, time and state units", {
  cp <- sqrt(2) - 1
  dp <- (1 + sqrt(5)) / 2
  for (k in 10^seq(-14, 14, 4)) {
    expect_equal(as.numeric(ctrl_care(-1, k, 1, k^2)$P), cp, tolerance = 1e-10)
    expect_equal(as.numeric(ctrl_dare(1, k, 1, k^2)$P), dp, tolerance = 1e-10)
    expect_equal(as.numeric(ctrl_care(-1, 1, k, k)$P), k * cp, tolerance = 1e-10)
    expect_equal(as.numeric(ctrl_care(-k, k, k, k)$P), cp, tolerance = 1e-10)
  }
  A <- matrix(c(-1, 1, -0.5, -2), 2, byrow = TRUE)
  B <- matrix(c(1, 0, 0.5, 1), 2, byrow = TRUE)
  Q <- matrix(c(2, 0.3, 0.3, 1), 2)
  R <- matrix(c(1, 0.2, 0.2, 2), 2)
  S <- matrix(c(0.1, 0.05, -0.05, 0.1), 2)
  base <- ctrl_care(A, B, Q, R, S)
  for (k in c(1e-12, 1e12)) {
    T <- diag(c(k, 1 / k))
    Ti <- diag(c(1 / k, k))
    x <- ctrl_care(Ti %*% A %*% T, Ti %*% B, T %*% Q %*% T, R, T %*% S)
    expect_equal(x$K, base$K %*% T, tolerance = 1e-10)
  }
})

test_that("Lyapunov factors match the dense solution", {
  set.seed(11)
  for (d in c(FALSE, TRUE)) {
    A <- matrix(rnorm(25), 5)
    A <- if (d) A / (1.2 * max(Mod(eigen(A)$values))) else
      A - (max(Re(eigen(A)$values)) + 0.5) * diag(5)
    B <- matrix(rnorm(10), 5)
    L <- controlcompR:::lyap_factor(A, B, d)
    X <- ctrl_lyap(A, B %*% t(B), discrete = d)
    expect_lt(max(abs(L %*% t(L) - X)) / max(abs(X)), 1e-12)
  }
})

test_that("well-conditioned poles of a stiff model are classified", {
  S <- ctrl_ss(diag(c(-1e-3, -1e11)), c(1e-3, 1e11), c(1, 1))
  expect_true(ctrl_is_stable(S))
  expect_equal(as.numeric(ctrl_norm(S, "inf")), 2, tolerance = 1e-8)
})

test_that("a tiny feedthrough gives a large finite zero", {
  expect_equal(Re(ctrl_zeros(ctrl_ss(-1, 1, 1, D = 1e-14))), -1e14 - 1,
               tolerance = 1e-10)
})

test_that("Hankel singular values of a model with complex poles match the Gramians", {
  set.seed(4)
  A <- -diag(c(0.5, 1, 2, 5, 10, 20, 50, 100)) + matrix(rnorm(64, sd = 0.3), 8)
  G <- ctrl_ss(A, matrix(rnorm(8), 8), matrix(rnorm(8), 1))
  hs <- ctrl_hsvd(G)
  ref <- sort(sqrt(Mod(eigen(ctrl_gram(G, "c") %*% ctrl_gram(G, "o"))$values)),
              decreasing = TRUE)
  expect_length(hs, 8)
  expect_equal(hs[1:4], ref[1:4], tolerance = 1e-8)
  expect_equal(attr(ctrl_balreal(G), "hsv")[1:4], ref[1:4], tolerance = 1e-8)
})

test_that("the reduction bound check rejects a reduced model that violates it", {
  G <- ctrl_ss(diag(c(-1, -2)), c(1, 1), c(1, 1))
  wrong <- ctrl_ss(-1, 1, 1)  # drops 1 / (s + 2), error 0.5 at DC
  expect_error(controlcompR:::check_reduction_bound(G, wrong, 1e-3, "truncate"),
               "certify")
  ok <- controlcompR:::check_reduction_bound(G, wrong, 1, "truncate")
  expect_identical(ok$status, "certified")
})

test_that("MHE matches the Kalman update for a rotated anisotropic prior", {
  U <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
  P0 <- U %*% diag(c(1, 1e-13)) %*% t(U)
  C <- matrix(c(0, 1), 1)
  sys <- ctrl_ss(diag(2), matrix(0, 2, 1), C, dt = 1)
  r <- ctrl_mhe_run(ctrl_mhe(sys, 1, Qn = P0, Rn = 1e-13, x0 = c(0, 0), P0 = P0),
                    1e-7)
  K <- P0 %*% t(C) %*% solve(C %*% P0 %*% t(C) + 1e-13)
  expect_equal(unname(r$estimates[1, ]), drop(K) * 1e-7, tolerance = 1e-6)
})

test_that("zeros from a tiny feedthrough keep the moderate zero accurate", {
  for (D in c(1e-7, 1e-10, 1e-14, 1e-30)) {
    z <- ctrl_zeros(ctrl_ss(diag(c(-1, -2)), c(1, 1), c(1, 1), D = D))
    ref <- polyroot(c(3 + 2 * D, 2 + 3 * D, D))
    ref <- ref[order(Mod(ref), decreasing = TRUE)]
    z <- z[order(Mod(z), decreasing = TRUE)]
    expect_equal(Re(z), Re(ref), tolerance = 1e-10)
  }
})

test_that("Riccati solutions are invariant to very different input units", {
  set.seed(2)
  A <- matrix(rnorm(9), 3)
  B <- matrix(rnorm(6), 3)
  S <- matrix(rnorm(6) * 0.1, 3)
  R <- diag(c(1, 2))
  Du <- diag(c(1e12, 1e-12))
  ref <- ctrl_care(A, B, diag(3), R, S)
  x <- ctrl_care(A, B %*% Du, diag(3), Du %*% R %*% Du, S %*% Du)
  expect_equal(x$P, ref$P, tolerance = 1e-12)
  expect_equal(Du %*% x$K, ref$K, tolerance = 1e-12)
  refd <- ctrl_dare(A / 3, B, diag(3), R, S)
  xd <- ctrl_dare(A / 3, B %*% Du, diag(3), Du %*% R %*% Du, S %*% Du)
  expect_equal(xd$P, refd$P, tolerance = 1e-12)
})

test_that("the UKF matches the Kalman filter for a rotated anisotropic prior", {
  U <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
  P0 <- U %*% diag(c(1, 1e-10)) %*% t(U)
  f <- ctrl_nlsys(function(k, x, u, p) x, function(k, x, u, p) sum(U[, 2] * x),
                  nstates = 2, ninputs = 0, dt = 1)
  ru <- ctrl_ukf(f, c(1e-5, 0), Qn = 0, Rn = 1e-12, x0 = c(0, 0), P0 = P0)
  kf <- ctrl_kalman_filter(ctrl_ss(diag(2), matrix(0, 2, 1), matrix(U[, 2], 1),
                                   dt = 1),
                           c(1e-5, 0), Qn = matrix(0, 2, 2), Rn = 1e-12, P0 = P0)
  expect_equal(unname(ru$filtered), unname(kf$filtered), tolerance = 1e-3)
})

test_that("frequency response of a badly scaled realization is evaluated", {
  for (a in c(1e100, 1e-100)) {
    s <- ctrl_tf2ss(ctrl_tf(1, c(1, a, 0, 0, 0)), minimal = FALSE)
    r <- drop(ctrl_freqresp(s, c(1, 0))$response)
    expect_equal(r[1], 1 / (1 - a * 1i), tolerance = 1e-12)
    expect_true(is.na(r[2]))  # w = 0 is a pole
    w <- 10^seq(-3, 3, length.out = 30)  # many points: Schur path
    r <- drop(ctrl_freqresp(s, w)$response)
    expect_equal(r, 1 / ((1i * w)^4 + a * (1i * w)^3), tolerance = 1e-12)
  }
})
