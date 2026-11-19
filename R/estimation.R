#' Luenberger observer
#'
#' Design an observer gain `L` by placing the eigenvalues of
#' \eqn{A - L C} with [ctrl_place()] applied to the dual pair
#' \eqn{(A^\top, C^\top)}. The observer is
#' \deqn{\dot{\hat x} = A \hat x + B u + L (y - C \hat x - D u)}
#' (or the corresponding difference equation in discrete time, a
#' one-step predictor).
#'
#' @param sys A `ctrl_ss` model (observable).
#' @param poles Desired observer poles.
#' @return An object of class `c("ctrl_observer", "ctrl_estimator")` with `L`,
#'   `poles` and `system` (the observer as a `ctrl_ss` model with inputs
#'   `u` then `y` and the state estimates as outputs).
#' @examples
#' plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
#' obs <- ctrl_observer(plant, c(-10, -11))
#' obs$L
#' @export
ctrl_observer <- function(sys, poles) {
  s <- as_ss(sys)
  if (!ctrl_is_observable(s)) {
    ctrl_abort("The model is not observable; observer poles cannot be ",
               "assigned freely.")
  }
  L <- t(ctrl_place(t(s$A), t(s$C), poles))
  structure(list(L = L, poles = attr(t(L), "poles") %||%
                   as.complex(eigen(s$A - L %*% s$C, only.values = TRUE)$values),
                 system = estimator_system(s, L), plant = s, dt = s$dt),
            class = c("ctrl_observer", "ctrl_estimator"))
}

estimator_system <- function(s, L) {
  n <- nrow(s$A)
  A <- s$A - L %*% s$C
  B <- cbind(s$B - L %*% s$D, L)
  ctrl_ss(A, B, diag(n), matrix(0, n, ncol(B)), dt = s$dt,
          input_names = c(s$input_names, s$output_names),
          output_names = paste0(s$state_names, "_hat"))
}

#' @export
print.ctrl_estimator <- function(x, digits = 4, ...) {
  cat(if (inherits(x, "ctrl_kalman")) "Steady-state Kalman filter" else
    "Luenberger observer", " (", time_domain_label(x$dt), ")\n", sep = "")
  cat("L =\n")
  print(signif(x$L, digits))
  cat("estimator poles: ", paste(format(signif(x$poles, digits)),
                                 collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Steady-state Kalman filter
#'
#' For the model
#' \deqn{x_{k+1} = A x_k + B u_k + G w_k, \quad y_k = C x_k + D u_k + v_k}
#' (or its continuous time analogue) with
#' \eqn{E[w w^\top] = Q_n}, \eqn{E[v v^\top] = R_n},
#' \eqn{E[w v^\top] = N_n}, compute the steady-state Kalman gains from the
#' filtering Riccati equation (the dual of the LQR equation).
#'
#' Discrete time conventions: `P` is the steady-state *prior* covariance
#' \eqn{P_{k|k-1}}; `M` is the measurement-update (innovation) gain
#' \eqn{M = P C^\top (C P C^\top + R_n)^{-1}} used in
#' \eqn{\hat x_{k|k} = \hat x_{k|k-1} + M (y_k - C\hat x_{k|k-1} - D u_k)};
#' `L` is the one-step predictor gain
#' \eqn{L = (A P C^\top + G N_n)(C P C^\top + R_n)^{-1}} used in
#' \eqn{\hat x_{k+1|k} = A \hat x_{k|k-1} + B u_k + L(y_k - C \hat x_{k|k-1}
#' - D u_k)}; `Z` is the posterior covariance \eqn{P_{k|k}} (for
#' \eqn{N_n = 0}). In continuous time, `L` is the Kalman-Bucy gain
#' \eqn{(P C^\top + G N_n) R_n^{-1}}.
#'
#' `Rn` must be positive definite: perfect (noise-free) measurements make
#' the filtering problem singular and are not silently regularized.
#' `Qn` may be singular (including zero) as long as the Riccati equation
#' has a stabilizing solution; this requires \eqn{(A, C)} detectable and no
#' unreachable noise modes on the stability boundary.
#'
#' @param sys A `ctrl_ss` model.
#' @param Qn Process noise covariance (size `ncol(G)`).
#' @param Rn Measurement noise covariance (`p x p`).
#' @param Nn Optional cross covariance (`ncol(G) x p`).
#' @param G Process noise input matrix; defaults to the identity (noise on
#'   every state).
#' @return An object of class `c("ctrl_kalman", "ctrl_estimator")` with `L`,
#'   `M` (discrete only), `P`, `Z` (discrete only), `poles` (eigenvalues of
#'   \eqn{A - L C}), and `system` (the predictor as a `ctrl_ss` model with
#'   inputs `u` then `y`).
#' @seealso [ctrl_kalman_filter()] for the time-varying filter.
#' @examples
#' sys <- ctrl_ss(matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE), c(0.005, 0.1),
#'                c(1, 0), dt = 0.1)
#' kf <- ctrl_kalman(sys, Qn = diag(c(1e-4, 1e-3)), Rn = 0.01)
#' kf$M
#' @export
ctrl_kalman <- function(sys, Qn, Rn, Nn = NULL, G = NULL) {
  s <- as_ss(sys)
  n <- nrow(s$A)
  p <- nrow(s$C)
  if (!n) ctrl_abort("Kalman filtering needs a model with states.")
  G <- if (is.null(G)) diag(n) else as_real_matrix(G, "G", nrow = n)
  check_dims(G, "G", n, ncol(G))
  q <- ncol(G)
  Qn <- check_psd(expand_weight(Qn, q, "Qn"), "Qn")
  Rn <- check_pd(expand_weight(Rn, p, "Rn"), "Rn")
  Nn <- if (is.null(Nn)) matrix(0, q, p) else as_real_matrix(Nn, "Nn", q, p)
  check_dims(Nn, "Nn", q, p)
  W0 <- rbind(cbind(Qn, Nn), cbind(t(Nn), Rn))
  if (min_eig_sym(W0) < -1e-10 * norm(W0, "2")) {
    ctrl_abort("The joint noise covariance [Qn Nn; Nn' Rn] is not positive ",
               "semidefinite.")
  }
  obs <- ctrl_observability(s)
  if (!obs$detectable) {
    ctrl_abort("(A, C) is not detectable; unobservable unstable modes: ",
               paste(format(signif(obs$unobservable_modes, 4)), collapse = ", "),
               ".")
  }
  Qx <- G %*% Qn %*% t(G)
  Sx <- G %*% Nn
  W <- rbind(cbind(Qx, Sx), cbind(t(Sx), Rn))
  if (min_eig_sym(W) < -1e-10 * norm(W, "2")) {
    ctrl_abort("The joint noise covariance [Qn Nn; Nn' Rn] is not positive ",
               "semidefinite.")
  }
  ric <- if (s$dt > 0) {
    ctrl_dare(t(s$A), t(s$C), Qx, Rn, Sx)
  } else {
    ctrl_care(t(s$A), t(s$C), Qx, Rn, Sx)
  }
  P <- ric$P
  if (s$dt > 0) {
    Sy <- s$C %*% P %*% t(s$C) + Rn
    M <- t(solve(Sy, s$C %*% P))
    L <- t(solve(Sy, s$C %*% P %*% t(s$A) + t(Sx)))
    Z <- symmetrize(P - M %*% Sy %*% t(M))
  } else {
    M <- NULL
    Z <- NULL
    L <- t(solve(Rn, s$C %*% P + t(Sx)))
  }
  dimnames(L) <- list(s$state_names, s$output_names)
  poles <- sort_poles(as.complex(eigen(s$A - L %*% s$C, only.values = TRUE)$values),
                      s$dt)
  structure(list(L = L, M = M, P = P, Z = Z, poles = poles, Qn = Qn, Rn = Rn,
                 Nn = Nn, G = G, system = estimator_system(s, L), plant = s,
                 dt = s$dt, residual = ric$residual),
            class = c("ctrl_kalman", "ctrl_estimator"))
}

#' Time-varying discrete Kalman filter and smoother
#'
#' Run the Kalman filter on measured data for
#' \eqn{x_{k+1} = A x_k + B u_k + G w_k}, \eqn{y_k = C x_k + D u_k + v_k}.
#'
#' **Conventions.** `x0` and `P0` are the prior mean and covariance of the
#' state at the first sample, \eqn{\hat x_{1|0}} and \eqn{P_{1|0}}. Each
#' step first *updates* with \eqn{y_k} and then *predicts* to \eqn{k + 1}:
#' \deqn{e_k = y_k - C \hat x_{k|k-1} - D u_k, \quad S_k = C P_{k|k-1}
#' C^\top + R_n, \quad K_k = P_{k|k-1} C^\top S_k^{-1},}
#' \deqn{\hat x_{k|k} = \hat x_{k|k-1} + K_k e_k, \quad
#' \hat x_{k+1|k} = A \hat x_{k|k} + B u_k.}
#' The covariance update uses the Joseph form
#' \eqn{P_{k|k} = (I - K_k C) P_{k|k-1} (I - K_k C)^\top + K_k R_n K_k^\top},
#' which keeps it symmetric positive semidefinite in floating point.
#' Rows of `y` containing `NA` are treated as missing measurements (the
#' update is skipped).
#'
#' `ctrl_kalman_smoother()` applies the Rauch-Tung-Striebel backward pass to
#' the filter output.
#'
#' @param sys A discrete time `ctrl_ss` model.
#' @param y Measurements, a matrix with one row per sample (`T x p`).
#' @param u Inputs (`T x m`), or `NULL` for zero input.
#' @param Qn,Rn Process and measurement noise covariances.
#' @param x0,P0 Prior mean and covariance at the first sample (defaults:
#'   zero and the identity).
#' @param G Process noise input matrix (default identity).
#' @return `ctrl_kalman_filter()` returns an object of class
#'   `ctrl_kalman_filter` with `filtered` (`T x n`), `predicted` (`T x n`,
#'   row `k` is \eqn{\hat x_{k|k-1}}), `P_filtered` and `P_predicted`
#'   (`n x n x T` arrays), `innovations`, `innovation_cov` (`p x p x T`),
#'   `gain` (`n x p x T`), `nis` (normalized innovation squared) and
#'   `loglik` (Gaussian log-likelihood of the measurements).
#'   `ctrl_kalman_smoother()` returns a list with `smoothed` and
#'   `P_smoothed`.
#' @examples
#' set.seed(1)
#' sys <- ctrl_ss(matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE), c(0.005, 0.1),
#'                c(1, 0), dt = 0.1)
#' x <- c(0, 1)
#' y <- numeric(50)
#' for (k in 1:50) {
#'   y[k] <- x[1] + rnorm(1, sd = 0.1)
#'   x <- drop(sys$A %*% x) + rnorm(2, sd = c(0.001, 0.01))
#' }
#' kf <- ctrl_kalman_filter(sys, y, Qn = diag(c(1e-6, 1e-4)), Rn = 0.01)
#' tail(kf$filtered, 3)
#' sm <- ctrl_kalman_smoother(kf)
#' @export
ctrl_kalman_filter <- function(sys, y, u = NULL, Qn, Rn, x0 = NULL, P0 = NULL,
                               G = NULL) {
  s <- as_ss(sys)
  if (s$dt == 0) {
    ctrl_abort("The time-varying Kalman filter needs a discrete time model; ",
               "discretize with ctrl_c2d() and the matching noise covariance.")
  }
  n <- nrow(s$A)
  p <- nrow(s$C)
  m <- ncol(s$B)
  Y <- if (is.null(dim(y))) matrix(as.numeric(y), ncol = 1) else as.matrix(y)
  if (ncol(Y) != p) ctrl_abort("`y` must have ", p, " columns.")
  if (any(is.infinite(Y)) || any(is.nan(Y))) {
    ctrl_abort("`y` contains infinite or NaN values; only NA marks missing data.")
  }
  Tn <- nrow(Y)
  U <- if (is.null(u)) matrix(0, Tn, m) else input_matrix(u, seq_len(Tn), m)
  G <- if (is.null(G)) diag(n) else as_real_matrix(G, "G", nrow = n)
  Qn <- check_psd(expand_weight(Qn, ncol(G), "Qn"), "Qn")
  Rn <- check_pd(expand_weight(Rn, p, "Rn"), "Rn")
  x <- if (is.null(x0)) rep(0, n) else as.numeric(x0)
  if (length(x) != n) ctrl_abort("`x0` must have length ", n, ".")
  P <- if (is.null(P0)) diag(n) else check_psd(expand_weight(P0, n, "P0"), "P0")
  Qx <- G %*% Qn %*% t(G)
  filt <- matrix(0, Tn, n)
  pred <- matrix(0, Tn, n)
  Pf <- array(0, c(n, n, Tn))
  Pp <- array(0, c(n, n, Tn))
  innov <- matrix(NA_real_, Tn, p)
  Sarr <- array(NA_real_, c(p, p, Tn))
  Karr <- array(0, c(n, p, Tn))
  nis <- rep(NA_real_, Tn)
  loglik <- 0
  In <- diag(n)
  for (k in seq_len(Tn)) {
    pred[k, ] <- x
    Pp[, , k] <- P
    yk <- Y[k, ]
    if (all(is.finite(yk))) {
      e <- yk - drop(s$C %*% x + s$D %*% U[k, ])
      S <- symmetrize(s$C %*% P %*% t(s$C) + Rn)
      cS <- chol(S)
      K <- t(backsolve(cS, forwardsolve(t(cS), s$C %*% P)))
      x <- x + drop(K %*% e)
      IKC <- In - K %*% s$C
      P <- symmetrize(IKC %*% P %*% t(IKC) + K %*% Rn %*% t(K))
      z <- forwardsolve(t(cS), e)
      nis[k] <- sum(z^2)
      loglik <- loglik - 0.5 * (p * log(2 * pi) + 2 * sum(log(diag(cS))) + nis[k])
      innov[k, ] <- e
      Sarr[, , k] <- S
      Karr[, , k] <- K
    } else if (any(is.finite(yk))) {
      ctrl_abort("Partially missing measurement rows are not supported (row ",
                 k, ").")
    }
    filt[k, ] <- x
    Pf[, , k] <- P
    x <- drop(s$A %*% x + s$B %*% U[k, ])
    P <- symmetrize(s$A %*% P %*% t(s$A) + Qx)
  }
  colnames(filt) <- colnames(pred) <- s$state_names
  structure(list(filtered = filt, predicted = pred, P_filtered = Pf,
                 P_predicted = Pp, innovations = innov, innovation_cov = Sarr,
                 gain = Karr, nis = nis, loglik = loglik, plant = s,
                 next_prediction = x, next_covariance = P, y = Y, u = U),
            class = "ctrl_kalman_filter")
}

#' @export
print.ctrl_kalman_filter <- function(x, ...) {
  cat("Kalman filter run over ", nrow(x$filtered), " samples, ",
      ncol(x$filtered), " states\n", sep = "")
  cat("  log-likelihood: ", signif(x$loglik, 6), "\n", sep = "")
  cat("  mean NIS: ", signif(mean(x$nis, na.rm = TRUE), 4), " (expected ",
      ncol(x$innovations), " for a consistent filter)\n", sep = "")
  invisible(x)
}

#' @rdname ctrl_kalman_filter
#' @param filter A `ctrl_kalman_filter` result.
#' @export
ctrl_kalman_smoother <- function(filter) {
  if (!inherits(filter, "ctrl_kalman_filter")) {
    ctrl_abort("`filter` must be a `ctrl_kalman_filter` result.")
  }
  if (inherits(filter, "ctrl_nl_filter")) {
    ctrl_abort("Smoothing is only implemented for the linear Kalman filter; ",
               "EKF and UKF results cannot be smoothed with ctrl_kalman_smoother().")
  }
  A <- filter$plant$A
  Tn <- nrow(filter$filtered)
  n <- ncol(filter$filtered)
  xs <- filter$filtered
  Ps <- filter$P_filtered
  if (Tn >= 2L) {
    for (k in (Tn - 1L):1) {
      Pf <- filter$P_filtered[, , k]
      Ppn <- filter$P_predicted[, , k + 1L]
      # pseudo-inverse: a singular predicted covariance (deterministic state
      # components) is valid and leaves those components unchanged
      J <- Pf %*% t(A) %*% MASS_ginv(Ppn, tol = 1e-12)
      xs[k, ] <- filter$filtered[k, ] + drop(J %*% (xs[k + 1L, ] -
                                                     filter$predicted[k + 1L, ]))
      Ps[, , k] <- symmetrize(Pf + J %*% (Ps[, , k + 1L] - Ppn) %*% t(J))
    }
  }
  list(smoothed = xs, P_smoothed = Ps)
}

#' @rdname ctrl_plots
#' @param truth Optional matrix of true states to overlay.
#' @export
plot.ctrl_kalman_filter <- function(x, truth = NULL, ...) {
  Tn <- nrow(x$filtered)
  rows <- list()
  for (j in seq_len(ncol(x$filtered))) {
    sd <- sqrt(pmax(x$P_filtered[j, j, ], 0))
    rows[[j]] <- data.frame(k = seq_len(Tn), state = colnames(x$filtered)[j],
                            estimate = x$filtered[, j],
                            lower = x$filtered[, j] - 2 * sd,
                            upper = x$filtered[, j] + 2 * sd,
                            truth = if (is.null(truth)) NA_real_ else
                              as.matrix(truth)[, j])
  }
  df <- do.call(rbind, rows)
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$k)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate, color = "estimate"))
  if (!is.null(truth)) {
    g <- g + ggplot2::geom_line(ggplot2::aes(y = .data$truth, color = "true"),
                                linetype = "dashed")
  }
  g + ggplot2::facet_wrap(~state, scales = "free_y") +
    ggplot2::labs(x = "sample", y = NULL, color = NULL,
                  title = "Kalman filter estimates (band: +/- 2 sd)") +
    ggplot2::theme_bw()
}

#' LQG controller
#'
#' Combine an LQR state feedback with a steady-state Kalman filter. Under
#' the separation principle (linear model, Gaussian noise with the given
#' covariances, and stabilizing LQR and Kalman designs) the closed-loop
#' poles are the union of the regulator poles, eig(A - BK), and the
#' estimator poles, eig(A - LC). LQG controllers have no guaranteed
#' robustness margins (Doyle, 1978); check margins of the actual loop.
#'
#' The controller maps measurements `y` to inputs `u` (it contains the
#' minus sign of the feedback):
#' * continuous time and discrete `form = "delayed"`:
#'   \eqn{\dot{\hat x} = (A - BK - LC + LDK)\hat x + L y}, \eqn{u = -K\hat x};
#' * discrete `form = "current"`: the input uses the measurement-updated
#'   estimate \eqn{\hat x_{k|k}}, \eqn{u_k = -K \hat x_{k|k}}, which reacts to
#'   \eqn{y_k} without a one-sample delay.
#'
#' @param sys Plant model.
#' @param Q,R,N LQR weights (see [ctrl_lqr()]).
#' @param Qn,Rn,Nn,G Noise model (see [ctrl_kalman()]).
#' @param form `"current"` or `"delayed"` (discrete models only).
#' @return An object of class `c("ctrl_lqg", "ctrl_controller")` with
#'   `controller` (a `ctrl_ss` model from `y` to `u`), `lqr`, `kalman` and
#'   `form`.
#' @references Doyle, J. C. (1978). Guaranteed margins for LQG regulators.
#'   *IEEE Transactions on Automatic Control*, 23(4), 756-757.
#' @examples
#' plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
#'                           c(1, 0)), 0.1)
#' lqg <- ctrl_lqg(plant, Q = diag(c(1, 0.1)), R = 0.1, Qn = diag(c(1e-4, 1e-3)),
#'                 Rn = 1e-2)
#' cl <- ctrl_closed_loop(plant, lqg)
#' ctrl_poles(cl)
#' @export
ctrl_lqg <- function(sys, Q, R, Qn, Rn, N = NULL, Nn = NULL, G = NULL,
                     form = c("current", "delayed")) {
  form <- match.arg(form)
  s <- as_ss(sys)
  lqr <- ctrl_lqr(s, Q, R, N)
  kf <- ctrl_kalman(s, Qn, Rn, Nn, G)
  K <- unname(lqr$K)
  L <- unname(kf$L)
  A <- s$A
  B <- s$B
  C <- s$C
  D <- s$D
  if (s$dt > 0 && form == "current") {
    M <- kf$M
    if (any(kf$Nn != 0)) {
      ctrl_abort("The current-estimator form assumes uncorrelated noise; use ",
                 "`form = \"delayed\"` with `Nn`.")
    }
    n <- nrow(A)
    W <- safe_solve(diag(ncol(B)) - K %*% M %*% D, diag(ncol(B)),
                    "current-estimator loop (I - K M D)")
    IMC <- diag(n) - M %*% C
    BA <- B - A %*% M %*% D
    Ac <- A %*% IMC - BA %*% W %*% K %*% IMC
    Bc <- A %*% M - BA %*% W %*% K %*% M
    Cc <- -W %*% K %*% IMC
    Dc <- -W %*% K %*% M
  } else {
    Ac <- A - B %*% K - L %*% C + L %*% D %*% K
    Bc <- L
    Cc <- -K
    Dc <- matrix(0, ncol(B), nrow(C))
  }
  ctrl <- ctrl_ss(Ac, Bc, Cc, Dc, dt = s$dt, input_names = s$output_names,
                  output_names = s$input_names,
                  state_names = paste0(s$state_names, "_hat"))
  structure(list(controller = ctrl, lqr = lqr, kalman = kf,
                 form = if (s$dt > 0) form else "continuous", plant = s,
                 dt = s$dt, K = lqr$K, L = kf$L),
            class = c("ctrl_lqg", "ctrl_controller"))
}

#' @export
print.ctrl_lqg <- function(x, digits = 4, ...) {
  cat("LQG controller (", time_domain_label(x$dt), ", ", x$form, " form)\n",
      sep = "")
  cat("regulator poles: ", paste(format(signif(x$lqr$poles, digits)),
                                 collapse = ", "), "\n", sep = "")
  cat("estimator poles: ", paste(format(signif(x$kalman$poles, digits)),
                                 collapse = ", "), "\n", sep = "")
  invisible(x)
}

lqg_closed_loop <- function(s, lqg, include_control) {
  ctrl <- lqg$controller
  cl <- ctrl_feedback(s, ctrl, sign = 1)
  Cy <- cl$C
  Dy <- cl$D
  n <- nrow(s$A)
  nc <- nrow(ctrl$A)
  out_names <- s$output_names
  C <- Cy
  D <- Dy
  if (include_control) {
    Cu <- cbind(matrix(0, ncol(s$B), n), ctrl$C) + ctrl$D %*% Cy
    Du <- diag(ncol(s$B)) + ctrl$D %*% Dy
    C <- rbind(C, Cu)
    D <- rbind(D, Du)
    out_names <- c(out_names, paste0("ctrl_", s$input_names))
  }
  ctrl_ss(cl$A, cl$B, C, D, dt = s$dt,
          state_names = c(s$state_names, ctrl$state_names),
          input_names = paste0("v_", s$input_names), output_names = out_names)
}

#' Stochastic closed-loop simulation of an LQG controller
#'
#' Simulate a discrete time plant with Gaussian process and measurement
#' noise under an LQG controller, optionally over several independent
#' realizations, and report the true and estimated state trajectories.
#' Summaries over realizations (mean squared estimation error, mean
#' quadratic cost) avoid judging a stochastic design from one favorable
#' trajectory.
#'
#' @param plant Discrete time plant (`ctrl_ss`); it may differ from the
#'   design model to study model mismatch, but must have the same
#'   dimensions.
#' @param lqg A `ctrl_lqg` object.
#' @param steps Number of samples.
#' @param x0 Initial true state.
#' @param Qn,Rn Covariances used to *generate* the noise (default: those of
#'   the design).
#' @param replications Number of independent realizations.
#' @param seed Optional seed (the previous RNG state is restored).
#' @return A list of class `ctrl_lqg_simulation` with arrays `x`, `xhat`,
#'   `u`, `y` (samples x signals x replications) and a `summary` data frame
#'   with the per-replication RMS estimation error and average stage cost.
#' @examples
#' plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
#'                           c(1, 0)), 0.1)
#' lqg <- ctrl_lqg(plant, Q = diag(c(1, 0.1)), R = 0.1, Qn = diag(c(1e-4, 1e-3)),
#'                 Rn = 1e-2)
#' sim <- ctrl_simulate_lqg(plant, lqg, steps = 100, x0 = c(1, 0),
#'                          replications = 5, seed = 42)
#' sim$summary
#' @export
ctrl_simulate_lqg <- function(plant, lqg, steps, x0, Qn = NULL, Rn = NULL,
                              replications = 1L, seed = NULL) {
  s <- as_ss(plant)
  if (!inherits(lqg, "ctrl_lqg")) ctrl_abort("`lqg` must be a `ctrl_lqg` object.")
  if (s$dt == 0) ctrl_abort("Stochastic simulation needs a discrete time plant.")
  if (!identical(dim(s$D), dim(lqg$plant$D)) || nrow(s$A) != nrow(lqg$plant$A)) {
    ctrl_abort("`plant` must have the dimensions of the design model.")
  }
  if (abs(s$dt - lqg$plant$dt) > 1e-12 * s$dt) {
    ctrl_abort("`plant` has sample interval ", s$dt, " but the controller was ",
               "designed for ", lqg$plant$dt, "; multirate loops are not supported.")
  }
  steps <- as_count(steps, "steps")
  replications <- as_count(replications, "replications")
  x0 <- as.numeric(x0)
  if (length(x0) != nrow(s$A) || any(!is.finite(x0))) {
    ctrl_abort("`x0` must be a finite vector of length ", nrow(s$A), ".")
  }
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  kf <- lqg$kalman
  Qn <- if (is.null(Qn)) kf$Qn else check_psd(expand_weight(Qn, ncol(kf$G), "Qn"), "Qn")
  Rn <- if (is.null(Rn)) kf$Rn else check_psd(expand_weight(Rn, nrow(s$C), "Rn"), "Rn")
  # process and measurement noise are drawn jointly with the design cross
  # covariance Nn
  q <- ncol(kf$G)
  Wj <- rbind(cbind(Qn, kf$Nn), cbind(t(kf$Nn), Rn))
  if (min_eig_sym(Wj) < -1e-10 * norm(Wj, "2")) {
    ctrl_abort("The joint noise covariance [Qn Nn; Nn' Rn] is not positive ",
               "semidefinite.")
  }
  Fj <- psd_sqrt_factor(Wj)
  n <- nrow(s$A)
  m <- ncol(s$B)
  p <- nrow(s$C)
  ctrl <- lqg$controller
  nc <- nrow(ctrl$A)
  X <- array(0, c(steps + 1L, n, replications))
  Xh <- array(0, c(steps + 1L, n, replications))
  U <- array(0, c(steps, m, replications))
  Yo <- array(0, c(steps, p, replications))
  Qw <- lqg$lqr$Q
  Rw <- lqg$lqr$R
  Nw <- lqg$lqr$N %||% matrix(0, n, m)
  summ <- data.frame(replication = seq_len(replications), rms_estimation_error = 0,
                     average_cost = 0)
  for (r in seq_len(replications)) {
    x <- as.numeric(x0)
    xc <- rep(0, nc)
    X[1, , r] <- x
    cost <- 0
    for (k in seq_len(steps)) {
      wv <- if (ncol(Fj)) drop(Fj %*% rnorm(ncol(Fj))) else rep(0, q + p)
      w <- wv[seq_len(q)]
      v <- wv[q + seq_len(p)]
      # u depends on y through Dc; with plant D the loop is algebraic
      E <- diag(p) - s$D %*% ctrl$D
      y <- drop(solve(E, s$C %*% x + s$D %*% (ctrl$C %*% xc) + v))
      u <- drop(ctrl$C %*% xc + ctrl$D %*% y)
      xhat <- if (lqg$form == "current") {
        Mk <- kf$M
        # the controller's own estimate uses the design measurement model
        xc + drop(Mk %*% (y - lqg$plant$C %*% xc - lqg$plant$D %*% u))
      } else xc
      Xh[k, , r] <- xhat
      U[k, , r] <- u
      Yo[k, , r] <- y
      cost <- cost + sum(x * (Qw %*% x)) + sum(u * (Rw %*% u)) +
        2 * sum(x * (Nw %*% u))
      x <- drop(s$A %*% x + s$B %*% u + kf$G %*% w)
      xc <- drop(ctrl$A %*% xc + ctrl$B %*% y)
      X[k + 1L, , r] <- x
    }
    Xh[steps + 1L, , r] <- if (lqg$form == "current") NA else xc
    err <- X[seq_len(steps), , r, drop = FALSE] - Xh[seq_len(steps), , r, drop = FALSE]
    summ$rms_estimation_error[r] <- sqrt(mean(err^2))
    summ$average_cost[r] <- cost / steps
  }
  structure(list(x = X, xhat = Xh, u = U, y = Yo, summary = summ,
                 time = (0:steps) * s$dt, state_names = s$state_names),
            class = "ctrl_lqg_simulation")
}
