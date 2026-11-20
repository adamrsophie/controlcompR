test_that("native Stein solver matches the reference block recursion", {
  set.seed(1)
  for (n in c(1, 2, 7, 30)) {
    A <- matrix(rnorm(n * n), n) / sqrt(n)
    A <- A / (max(Mod(eigen(A)$values)) * 1.1)
    sch <- controlcompR:::schur_real(A)
    Q <- crossprod(matrix(rnorm(n * n), n))
    Qt <- t(sch$Z) %*% Q %*% sch$Z
    expect_equal(controlcompR:::stein_schur(sch$T, Qt),
                 controlcompR:::stein_schur_r(sch$T, Qt), tolerance = 1e-12)
  }
})

test_that("Schur-based and dense frequency responses agree", {
  set.seed(2)
  n <- 12
  G <- ctrl_ss(matrix(rnorm(n * n), n) / 3 - diag(n), matrix(rnorm(2 * n), n),
               matrix(rnorm(2 * n), 2), matrix(rnorm(4), 2))
  w <- 10^seq(-2, 2, length.out = 40)
  many <- ctrl_freqresp(G, w)$response
  few <- vapply(w[c(1, 20, 40)], function(wi) ctrl_freqresp(G, wi)$response,
                array(0i, c(2, 2, 1)))
  expect_equal(many[, , c(1, 20, 40)], array(few, c(2, 2, 3)), tolerance = 1e-12)
})

test_that("filters can keep only the variances", {
  sys <- ctrl_ss(matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE), c(0, 0.1), c(1, 0), dt = 0.1)
  set.seed(3)
  y <- cumsum(rnorm(30, sd = 0.1))
  full <- ctrl_kalman_filter(sys, y, Qn = diag(2) * 1e-3, Rn = 0.01)
  lean <- ctrl_kalman_filter(sys, y, Qn = diag(2) * 1e-3, Rn = 0.01,
                             covariances = "variances")
  expect_null(lean$P_filtered)
  expect_equal(lean$filtered, full$filtered)
  expect_equal(unname(lean$variances[, 2]), full$P_filtered[2, 2, ])
  expect_error(ctrl_kalman_smoother(lean), "covariance history")
  f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2])
  ns <- ctrl_nlsys(f, function(k, x, u, p) x[1], nstates = 2, ninputs = 0, dt = 0.1)
  e1 <- ctrl_ekf(ns, y, Qn = diag(2) * 1e-3, Rn = 0.01)
  e2 <- ctrl_ukf(ns, y, Qn = diag(2) * 1e-3, Rn = 0.01, covariances = "variances")
  expect_null(e2$P_filtered)
  expect_equal(unname(e1$variances[, 1]), e1$P_filtered[1, 1, ])
})

test_that("analytic model Jacobians are validated and used", {
  f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] + 0.1 * (-sin(x[1]) + u))
  J <- function(k, x, u, p) {
    list(A = matrix(c(1, -0.1 * cos(x[1]), 0.1, 1), 2), B = matrix(c(0, 0.1), 2))
  }
  expect_error(ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1,
                          jacobian = function(...) list(A = diag(3))), "2 x 2")
  expect_error(ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1,
                          jacobian = function(...) list(Z = 1)), "named list")
  sj <- ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1, jacobian = J)
  lin <- ctrl_linearize(sj, c(0.3, 0), 0)
  expect_identical(attr(lin, "method"), "analytic (user supplied)")
  expect_equal(lin$A, J(0, c(0.3, 0), 0, NULL)$A, ignore_attr = TRUE)
  sf <- ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1)
  a <- ctrl_nmpc_step(ctrl_nmpc(sj, 10, diag(c(10, 1)), 0.1, u_min = -2, u_max = 2), c(1, 0))
  b <- ctrl_nmpc_step(ctrl_nmpc(sf, 10, diag(c(10, 1)), 0.1, u_min = -2, u_max = 2), c(1, 0))
  expect_equal(a$u_plan, b$u_plan, tolerance = 1e-6)
})

test_that("ARX order selection with a shared regressor bank matches single fits", {
  set.seed(4)
  u <- rnorm(300)
  y <- as.numeric(stats::filter(0.5 * c(0, u[-300]), 0.8, "recursive")) + rnorm(300, sd = 0.05)
  d <- ctrl_iddata(y, u)
  tab <- ctrl_arx_order(d, 1:3, 1:3)
  start <- 4L
  for (r in seq_len(nrow(tab))) {
    m <- ctrl_arx(d, tab$na[r], tab$nb[r], tab$nk[r], start = start)
    expect_equal(tab$bic[r], m$bic)
  }
})

test_that("fixed-lag smoothing, final storage and chunked filtering agree", {
  set.seed(3)
  sys <- ctrl_ss(matrix(c(1, 0.1, 0, 0.95), 2, byrow = TRUE), c(0.005, 0.1),
                 c(1, 0), dt = 0.1)
  Tn <- 40
  y <- cumsum(rnorm(Tn, sd = 0.1))
  y[c(10, 11)] <- NA
  Qn <- diag(c(1e-4, 1e-3))
  kf <- function(yy, ...) ctrl_kalman_filter(sys, yy, Qn = Qn, Rn = 0.01, ...)
  full <- kf(y)
  expect_equal(kf(y, lag = Tn)$lagged, ctrl_kalman_smoother(full)$smoothed)
  expect_equal(kf(y, lag = 0)$lagged, full$filtered)
  L <- 4
  fl <- kf(y, lag = L, covariances = "final")
  for (k in c(1, 9, 20, Tn - 2, Tn)) {
    ref <- ctrl_kalman_smoother(kf(y[seq_len(min(k + L, Tn))]))$smoothed[k, ]
    expect_equal(fl$lagged[k, ], ref)
  }
  expect_null(fl$variances)
  expect_null(fl$P_filtered)
  expect_equal(fl$P_final, full$P_filtered[, , Tn])
  expect_error(kf(y, lag = -1), "nonnegative")
  expect_error(ctrl_kalman_smoother(fl), "covariance history")
  a <- kf(y[1:25])
  b <- kf(y[26:Tn], x0 = a$next_prediction, P0 = a$next_covariance)
  expect_equal(rbind(a$filtered, b$filtered), full$filtered)
  expect_equal(a$loglik + b$loglik, full$loglik)
  nl <- ctrl_nlsys(update = function(t, x, u, p) 0.9 * x,
                   output = function(t, x, u, p) x, nstates = 1, ninputs = 0, dt = 1)
  e <- ctrl_ekf(nl, y, Qn = 0.01, Rn = 0.01, covariances = "final")
  expect_null(e$variances)
  expect_equal(drop(e$P_final), ctrl_ekf(nl, y, Qn = 0.01, Rn = 0.01)$P_filtered[1, 1, Tn])
})

test_that("Householder staircase keeps the controllable subspace and similarity", {
  set.seed(7)
  n <- 30
  A <- matrix(rnorm(n * n), n) / sqrt(n)
  B <- matrix(rnorm(n), n)
  # make (A, B) uncontrollable: a 10 state block not reached by B
  A[21:30, 1:20] <- 0
  B[21:30, ] <- 0
  st <- staircase(A, B)
  expect_identical(st$nc, 20L)
  expect_equal(crossprod(st$U), diag(n))
  expect_equal(t(st$U) %*% A %*% st$U, st$A)
  expect_equal(max(abs(st$A[21:30, 1:20])), 0, tolerance = 1e-12)
  expect_equal(max(abs(st$B[21:30, ])), 0, tolerance = 1e-12)
})

test_that("N4SID with the native QR and direct Hankel stack is unchanged", {
  set.seed(2)
  true <- ctrl_ss(matrix(c(0.7, 0.2, -0.2, 0.7), 2), c(1, 0.5), c(1, 1), dt = 1)
  u <- rnorm(600)
  y <- ctrl_simulate(true, input = u, times = 0:599)$y[, 1, 1]
  m <- ctrl_n4sid(ctrl_iddata(y, u), order = 2, feedthrough = TRUE)
  expect_equal(sort(Mod(ctrl_poles(m$sys))), sort(Mod(ctrl_poles(true))), tolerance = 1e-8)
  expect_lt(max(abs(m$residuals)), 1e-8)
  expect_error(ctrl_n4sid(ctrl_iddata(y, rep(1, 600)), order = 2), "persistently exciting")
})
