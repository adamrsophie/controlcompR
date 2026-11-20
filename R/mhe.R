#' Moving horizon estimation for linear models
#'
#' Estimate the state of
#' \eqn{x_{k+1} = A x_k + B u_k + G w_k}, \eqn{y_k = C x_k + D u_k + v_k}
#' by solving, at each sample \eqn{T}, the finite-horizon problem over the
#' window \eqn{s = T - N + 1, \ldots, T}:
#' \deqn{\min_{x_s, \ldots, x_T,\, w_s, \ldots, w_{T-1}}
#' \|x_s - \bar x_s\|^2_{P_s^{-1}} + \sum_{j=s}^{T-1} \|w_j\|^2_{Q_n^{-1}} +
#' \sum_{j=s}^{T} \|y_j - C x_j - D u_j\|^2_{R_n^{-1}}}
#' subject to the model equations and optional bounds on the states and the
#' process noise. The estimate is \eqn{\hat x_{T|T} = x_T}.
#'
#' **Arrival cost.** The first term summarizes the data before the window.
#' The prior mean is the prediction from the previous MHE estimate,
#' \eqn{\bar x_s = A \hat x_{s-1|s-1} + B u_{s-1}}, and \eqn{P_s} follows the
#' Kalman filter covariance recursion. Without active constraints this
#' makes MHE coincide with the Kalman filter (which the tests verify); with
#' active constraints MHE differs from, and can improve on, the Kalman
#' filter, at the price of one QP per sample, so its cost grows with the
#' horizon. MHE is worth it when states or noises have hard physical bounds
#' (concentrations, saturations) that the Kalman filter would violate.
#'
#' Each window QP is solved by OSQP. Once the window is full, the solver
#' workspace is reused from sample to sample while the sparsity pattern of
#' the QP does not change; only its numbers are updated.
#'
#' @param sys Discrete time `ctrl_ss` model.
#' @param horizon Window length `N` (number of measurements in the window).
#' @param Qn,Rn Process and measurement noise covariances (positive
#'   definite).
#' @param G Process noise input matrix (default identity).
#' @param x0,P0 Prior mean and covariance at the first sample.
#' @param x_min,x_max Bounds on the states.
#' @param w_min,w_max Bounds on the process noise.
#' @param solver_settings OSQP settings; they override the `accuracy`
#'   profile.
#' @param accuracy Solver accuracy profile, as in [ctrl_mpc()] (`"high"`
#'   uses tolerances 1e-10).
#' @return An object of class `ctrl_mhe`; run it with `ctrl_mhe_run()`.
#' @examples
#' sys <- ctrl_ss(0.9, 1, 1, dt = 1)
#' mhe <- ctrl_mhe(sys, horizon = 5, Qn = 0.01, Rn = 0.1, x_min = 0)
#' est <- ctrl_mhe_run(mhe, y = c(0.2, -0.1, 0.4, 0.3, 0.1, 0.5))
#' est$estimates
#' @export
ctrl_mhe <- function(sys, horizon, Qn, Rn, G = NULL, x0 = NULL, P0 = NULL,
                     x_min = -Inf, x_max = Inf, w_min = -Inf, w_max = Inf,
                     solver_settings = list(),
                     accuracy = c("high", "standard", "fast")) {
  accuracy <- match.arg(accuracy)
  s <- as_ss(sys)
  if (s$dt == 0) ctrl_abort("MHE needs a discrete time model.")
  n <- nrow(s$A)
  p <- nrow(s$C)
  N <- as_count(horizon, "horizon")
  G <- if (is.null(G)) diag(n) else as_real_matrix(G, "G", nrow = n)
  q <- ncol(G)
  Qn <- check_pd(expand_weight(Qn, q, "Qn"), "Qn")
  Rn <- check_pd(expand_weight(Rn, p, "Rn"), "Rn")
  x0 <- if (is.null(x0)) rep(0, n) else as.numeric(x0)
  P0 <- if (is.null(P0)) diag(n) else check_pd(expand_weight(P0, n, "P0"), "P0")
  structure(list(plant = s, horizon = N, Qn = Qn, Rn = Rn, G = G, x0 = x0,
                 P0 = P0, x_min = expand_bound(x_min, n, "x_min"),
                 x_max = expand_bound(x_max, n, "x_max"),
                 w_min = expand_bound(w_min, q, "w_min"),
                 w_max = expand_bound(w_max, q, "w_max"),
                 settings = utils::modifyList(osqp_profile(accuracy, 1e-10),
                                              solver_settings)),
            class = "ctrl_mhe")
}

#' Run a moving horizon estimator on data
#'
#' @param mhe A `ctrl_mhe` object.
#' @param y Measurements (`T x p`, `NA` rows are missing samples).
#' @param u Inputs (`T x m`) or `NULL`.
#' @param warm_start Start each window's QP from the previous window's
#'   solution, shifted by one sample (the new last state is the model
#'   prediction and the new last disturbance is zero). This changes only
#'   the number of OSQP iterations, not the solution.
#' @return A list of class `ctrl_mhe_result` with `estimates` (`T x n`,
#'   row `k` is \eqn{\hat x_{k|k}}), `status` (OSQP status per sample),
#'   `objective`, `iterations` and `solve_time`.
#' @rdname ctrl_mhe
#' @export
ctrl_mhe_run <- function(mhe, y, u = NULL, warm_start = TRUE) {
  if (!inherits(mhe, "ctrl_mhe")) ctrl_abort("`mhe` must be a `ctrl_mhe` object.")
  s <- mhe$plant
  n <- nrow(s$A)
  m <- ncol(s$B)
  p <- nrow(s$C)
  q <- ncol(mhe$G)
  Y <- measurement_matrix(y, p)
  Tn <- nrow(Y)
  U <- if (is.null(u)) matrix(0, Tn, m) else input_matrix(u, seq_len(Tn), m)
  Ri <- solve(mhe$Rn)
  Qi <- solve(mhe$Qn)
  Qx <- mhe$G %*% mhe$Qn %*% t(mhe$G)
  est <- matrix(NA_real_, Tn, n, dimnames = list(NULL, s$state_names))
  status <- character(Tn)
  obj <- numeric(Tn)
  stime <- numeric(Tn)
  iters <- integer(Tn)
  zprev <- NULL
  st_prev <- 0L
  model <- Pprev <- Aprev <- NULL
  # prior covariance at the window start: the Kalman prediction recursion,
  # advanced lazily (the window start never moves back), so only the
  # current covariance is kept instead of one n x n matrix per sample
  P <- mhe$P0
  p_at <- 1L
  advance_prior <- function(to) {
    while (p_at < to) {
      if (all(is.finite(Y[p_at, ]))) {
        S <- s$C %*% P %*% t(s$C) + mhe$Rn
        K <- P %*% t(s$C) %*% solve(S)
        IKC <- diag(n) - K %*% s$C
        P <<- symmetrize(IKC %*% P %*% t(IKC) + K %*% mhe$Rn %*% t(K))
      }
      P <<- symmetrize(s$A %*% P %*% t(s$A) + Qx)
      p_at <<- p_at + 1L
    }
    P
  }
  for (Tk in seq_len(Tn)) {
    st <- max(1L, Tk - mhe$horizon + 1L)
    L <- Tk - st + 1L
    xbar <- if (st == 1L) mhe$x0 else
      drop(s$A %*% est[st - 1L, ] + s$B %*% U[st - 1L, ])
    # arrival cost on the support of the prior covariance; components with
    # zero prior variance are fixed to the prior mean. The support is
    # decided on the equilibrated covariance, so a small positive variance
    # stays a (strong) prior and is not turned into an equality.
    sp <- psd_pinv(advance_prior(st))
    Pinv <- sp$inv
    Nul <- sp$null
    const <- sum(xbar * (Pinv %*% xbar))
    xi <- function(j) (j - 1L) * n + seq_len(n)
    wi <- function(j) L * n + (j - 1L) * q + seq_len(q)
    nz <- L * n + (L - 1L) * q
    Hm <- trip_new()
    qv <- numeric(nz)
    trip_add(Hm, xi(1), xi(1), 2 * Pinv)
    qv[xi(1)] <- -2 * drop(Pinv %*% xbar)
    for (j in seq_len(L)) {
      k <- st + j - 1L
      if (all(is.finite(Y[k, ]))) {
        trip_add(Hm, xi(j), xi(j), 2 * t(s$C) %*% Ri %*% s$C)
        rk <- drop(Y[k, ] - s$D %*% U[k, ])
        qv[xi(j)] <- qv[xi(j)] - 2 * drop(t(s$C) %*% Ri %*% rk)
        const <- const + sum(rk * (Ri %*% rk))
      }
      if (j < L) trip_add(Hm, wi(j), wi(j), 2 * Qi)
    }
    Ac <- trip_new()
    lo <- hi <- list()
    r <- 0L
    if (L > 1L) {
      for (j in seq_len(L - 1L)) {
        rows <- r + seq_len(n)
        trip_add(Ac, rows, xi(j + 1L), diag(n))
        trip_add(Ac, rows, xi(j), -s$A)
        trip_add(Ac, rows, wi(j), -mhe$G)
        bu <- drop(s$B %*% U[st + j - 1L, ])
        lo[[length(lo) + 1L]] <- bu
        hi[[length(hi) + 1L]] <- bu
        r <- r + n
      }
    }
    if (ncol(Nul)) {
      rows <- r + seq_len(ncol(Nul))
      trip_add(Ac, rows, xi(1), t(Nul))
      nb <- drop(t(Nul) %*% xbar)
      lo[[length(lo) + 1L]] <- nb
      hi[[length(hi) + 1L]] <- nb
      r <- r + ncol(Nul)
    }
    if (any(is.finite(c(mhe$x_min, mhe$x_max)))) {
      for (j in seq_len(L)) {
        rows <- r + seq_len(n)
        trip_add(Ac, rows, xi(j), diag(n))
        lo[[length(lo) + 1L]] <- rep_len(mhe$x_min, n)
        hi[[length(hi) + 1L]] <- rep_len(mhe$x_max, n)
        r <- r + n
      }
    }
    if (L > 1L && any(is.finite(c(mhe$w_min, mhe$w_max)))) {
      for (j in seq_len(L - 1L)) {
        rows <- r + seq_len(q)
        trip_add(Ac, rows, wi(j), diag(q))
        lo[[length(lo) + 1L]] <- rep_len(mhe$w_min, q)
        hi[[length(hi) + 1L]] <- rep_len(mhe$w_max, q)
        r <- r + q
      }
    }
    if (r == 0L) {
      # OSQP needs at least one constraint row; add a free one
      trip_add(Ac, 1L, 1L, matrix(1))
      lo <- list(-Inf)
      hi <- list(Inf)
      r <- 1L
    }
    lo <- unlist(lo, use.names = FALSE)
    hi <- unlist(hi, use.names = FALSE)
    H <- trip_sparse(Hm, c(nz, nz))
    A <- trip_sparse(Ac, c(r, nz))
    # change of variables z = M zs before calling the solver, so that
    # solver tolerances do not depend on the units or the orientation of
    # the prior: the first state of the window is written in the prior's
    # own coordinates [F N] (F F' = P on its support, N its deterministic
    # directions), which turns the arrival cost into an identity block even
    # for a rotated, very anisotropic covariance; then every variable is
    # scaled to unit curvature and the objective is normalized. The
    # solution, objective and duals are mapped back afterwards (constraint
    # rows are not rescaled, so the duals keep their meaning).
    T1 <- cbind(sp$factor, Nul)
    Mt <- Matrix::bdiag(Matrix::Matrix(T1, sparse = TRUE), Matrix::Diagonal(nz - n))
    H1 <- Matrix::t(Mt) %*% H %*% Mt
    hd <- Matrix::diag(H1)
    dv <- ifelse(is.finite(hd) & hd > 0, 1 / sqrt(hd), 1)
    Mv <- Mt %*% Matrix::Diagonal(x = dv)
    Hs <- Matrix::t(Mv) %*% H %*% Mv
    Hs <- (Hs + Matrix::t(Hs)) / 2
    osc <- max(abs(Matrix::diag(Hs)))
    if (!is.finite(osc) || osc <= 0) osc <- 1
    Pm <- methods::as(Matrix::triu(Hs / osc), "CsparseMatrix")
    Am <- methods::as(A %*% Mv, "CsparseMatrix")
    qs <- as.numeric(Matrix::t(Mv) %*% qv) / osc
    # once the window is full the QP keeps its sparsity pattern from sample
    # to sample (unless missing data or a change of support alters it), so
    # the solver workspace is reused and only the numbers are updated
    if (!is.null(model) && same_pattern(Pm, Pprev) && same_pattern(Am, Aprev)) {
      model@Update(q = qs, l = lo, u = hi, Px = Pm@x, Ax = Am@x)
      if (!isTRUE(warm_start)) model@ColdStart()
    } else {
      model <- osqp::osqp(P = Pm, q = qs, A = Am, l = lo, u = hi,
                          pars = do.call(osqp::osqpSettings, mhe$settings))
    }
    Pprev <- Pm
    Aprev <- Am
    to_scaled <- function(z) {
      z[seq_len(n)] <- solve(T1, z[seq_len(n)])
      z / dv
    }
    lay <- c(dyn = L - 1L, nul = ncol(Nul),
             xb = if (any(is.finite(c(mhe$x_min, mhe$x_max)))) L else 0L,
             wb = if (L > 1L && any(is.finite(c(mhe$w_min, mhe$w_max)))) L - 1L else 0L,
             free = as.integer(r == 1L && L == 1L && !ncol(Nul) &&
                                 !any(is.finite(c(mhe$x_min, mhe$x_max)))))
    if (isTRUE(warm_start) && !is.null(zprev)) {
      sh <- st - st_prev
      yw <- c(shift_blocks(yprev, n, lay_prev[["dyn"]], lay[["dyn"]], sh),
              rep(0, lay[["nul"]]),
              shift_blocks(drop_head(yprev, n * lay_prev[["dyn"]] + lay_prev[["nul"]]),
                           n, lay_prev[["xb"]], lay[["xb"]], sh),
              shift_blocks(drop_head(yprev, n * (lay_prev[["dyn"]] + lay_prev[["xb"]]) +
                                       lay_prev[["nul"]]),
                           q, lay_prev[["wb"]], lay[["wb"]], sh),
              rep(0, lay[["free"]]))
      model@WarmStart(x = to_scaled(mhe_shift(zprev, sh, L, n, q, s, U[Tk - 1L, ])),
                      y = if (length(yw) == r) yw / osc else NULL)
    }
    res <- osqp_solve(model)
    status[Tk] <- res$info$status
    if (!identical(res$info$status, "solved")) {
      ctrl_abort("MHE QP failed at sample ", Tk, " with status \"",
                 res$info$status, "\".")
    }
    zsol <- as.numeric(Mv %*% res$x)
    est[Tk, ] <- zsol[xi(L)]
    obj[Tk] <- osc * res$info$obj_val + const
    stime[Tk] <- res$info$solve_time
    iters[Tk] <- as.integer(res$info$iter %||% NA_integer_)
    zprev <- zsol
    yprev <- osc * res$y
    lay_prev <- lay
    st_prev <- st
  }
  structure(list(estimates = est, status = status, objective = obj,
                 iterations = iters, solve_time = stime),
            class = "ctrl_mhe_result")
}

# Warm start for the next window from the previous solution z (states
# x_1..x_Lp, then disturbances w_1..w_{Lp-1}). The window start moved by
# `shift` (0 while the window grows, 1 once it is full); the new last state
# is the model prediction and the new last disturbance is zero.
mhe_shift <- function(z, shift, L, n, q, s, u_last) {
  Lp <- L - 1L + shift
  X <- matrix(z[seq_len(Lp * n)], n)
  W <- if (Lp > 1L) matrix(z[Lp * n + seq_len((Lp - 1L) * q)], q) else matrix(0, q, 0)
  xn <- drop(s$A %*% X[, Lp] + s$B %*% u_last)
  X <- cbind(X, xn)[, shift + seq_len(L), drop = FALSE]
  W <- cbind(W, rep(0, q))[, shift + seq_len(L - 1L), drop = FALSE]
  c(as.vector(X), as.vector(W))
}

# Blocks of `size` entries: drop the first `shift`, keep `new` blocks,
# padding with zeros.
shift_blocks <- function(v, size, prev, new, shift) {
  if (!new) return(numeric())
  b <- if (prev) v[seq_len(prev * size)] else numeric()
  b <- c(drop_head(b, min(shift, prev) * size), rep(0, new * size))
  b[seq_len(new * size)]
}

# Drop the first k entries; v[-seq_len(0)] would drop everything.
drop_head <- function(v, k) if (k > 0) v[-seq_len(k)] else v

# TRUE when two sparse (CSC) matrices have the same dimensions and pattern.
same_pattern <- function(X, Y) {
  !is.null(Y) && identical(X@Dim, Y@Dim) && identical(X@p, Y@p) &&
    identical(X@i, Y@i)
}
