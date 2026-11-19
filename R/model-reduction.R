#' Minimal realization
#'
#' Remove uncontrollable and unobservable states with two orthogonal
#' staircase reductions (see [ctrl_controllability()]). The input-output
#' behavior is preserved up to the rank tolerance.
#'
#' @param sys A model.
#' The default tolerance first tries
#' \eqn{\sqrt{\epsilon}\max(\|A\|, \|B\|, \|C\|)} (as in MATLAB), which
#' removes modes that cancel to roughly half the working precision, as
#' needed for realizations built from polynomial coefficients. The reduced
#' model is then compared with the original at frequencies spread over the
#' range of the poles; when they disagree (a badly scaled model whose small
#' couplings were mistaken for zeros), the reduction is repeated with the
#' much tighter staircase tolerance \eqn{10 n^2 \epsilon \|[A, B]\|}, and
#' an error is raised if even that result does not match.
#'
#' @param tol Absolute rank tolerance for the staircase reductions. When
#'   given, it is used as is and the frequency check is skipped.
#' @return A `ctrl_ss` model with a minimal number of states.
#' @examples
#' # Two copies of the same first-order system: one state is redundant
#' sys <- ctrl_ss(diag(c(-1, -1)), c(1, 1), c(1, 1))
#' ctrl_dims(ctrl_minreal(sys))
#' @export
ctrl_minreal <- function(sys, tol = NULL) {
  if (is_ctrl_tf(sys) && length(sys$num) == 1L) {
    sys <- ctrl_tf2ss(sys, minimal = FALSE)
  }
  s <- as_ss(sys)
  n <- nrow(s$A)
  if (!n) return(s)
  if (!is.null(tol)) return(minreal_tol(s, tol))
  loose <- sqrt(.Machine$double.eps) *
    max(norm(s$A, "2"), norm(s$B, "2"), norm(s$C, "2"), 1)
  red <- minreal_tol(s, loose)
  if (same_response(s, red)) return(red)
  tight <- 10 * n^2 * .Machine$double.eps *
    max(norm(s$A, "2"), norm(s$B, "2"), norm(s$C, "2"), 1)
  red <- minreal_tol(s, tight)
  if (same_response(s, red)) return(red)
  ctrl_abort("The minimal realization does not reproduce the frequency ",
             "response of the model; the model is too badly scaled for a ",
             "reliable reduction. Rescale states, inputs or outputs.")
}

# Compare two realizations at frequencies spread over the pole magnitudes.
same_response <- function(s, r, rtol = 1e-6) {
  if (nrow(r$A) == nrow(s$A)) return(TRUE)
  ev <- Mod(eigen(s$A, only.values = TRUE)$values)
  if (s$dt > 0) {
    z <- c(1, -1, exp(1i * c(0.001, 0.01, 0.1, 0.5, 1, 2, 3) * 1.37))
  } else {
    ev <- ev[ev > 0]
    lo <- if (length(ev)) min(ev) else 1
    hi <- if (length(ev)) max(ev) else 1
    w <- 10^seq(log10(lo) - 2, log10(hi) + 2, length.out = 25) * 1.0137
    z <- 1i * c(0, w)
  }
  g1 <- eval_system(s, z)
  g2 <- eval_system(r, z)
  ok <- is.finite(Re(g1)) & is.finite(Re(g2))
  if (!any(ok)) return(TRUE)
  scale <- max(Mod(g1[ok]), Mod(s$D), .Machine$double.xmin)
  # compare per frequency relative to the response size there and overall
  err <- Mod(g1 - g2)[ok]
  ref <- pmax(Mod(g1)[ok], 1e-8 * scale)
  all(err <= rtol * ref + 1e-10 * scale)
}

minreal_tol <- function(s, tol) {
  st <- staircase(s$A, s$B, tol)
  nc <- st$nc
  keep <- seq_len(nc)
  Ac <- st$A[keep, keep, drop = FALSE]
  Bc <- st$B[keep, , drop = FALSE]
  Cc <- s$C %*% st$U[, keep, drop = FALSE]
  if (nc == 0L) {
    return(ctrl_gain(s$D, dt = s$dt, input_names = s$input_names,
                     output_names = s$output_names))
  }
  so <- staircase(t(Ac), t(Cc), tol)
  no <- so$nc
  if (no == 0L) {
    return(ctrl_gain(s$D, dt = s$dt, input_names = s$input_names,
                     output_names = s$output_names))
  }
  ko <- seq_len(no)
  Ao <- t(so$A)[ko, ko, drop = FALSE]
  Co <- t(so$B)[, ko, drop = FALSE]
  Bo <- t(so$U)[ko, , drop = FALSE] %*% Bc
  ctrl_ss(Ao, Bo, Co, s$D, dt = s$dt, input_names = s$input_names,
          output_names = s$output_names)
}

#' Stable and unstable decomposition
#'
#' Split a model into `G = G_stable + G_unstable` with an ordered real Schur
#' form followed by a Sylvester-equation block diagonalization. The feed
#' through `D` is assigned to the stable part. Poles on the stability
#' boundary (within `tol`) make the split ill-defined and cause an error.
#'
#' @param sys A model.
#' @param tol Boundary tolerance.
#' @return A list with `stable` and `unstable` `ctrl_ss` models.
#' @examples
#' sys <- ctrl_ss(diag(c(-1, 2, -3)), c(1, 1, 1), c(1, 1, 1))
#' parts <- ctrl_stabsep(sys)
#' ctrl_poles(parts$unstable)
#' @export
ctrl_stabsep <- function(sys, tol = NULL) {
  s <- as_ss(sys)
  n <- nrow(s$A)
  st <- ctrl_stability(s, tol = tol)
  if (length(st$boundary_poles)) {
    ctrl_abort("The model has poles on the stability boundary; the stable / ",
               "unstable decomposition is not defined.")
  }
  if (!n) return(list(stable = s, unstable = ctrl_gain(0 * s$D, dt = s$dt)))
  sch <- schur_real(s$A, order = if (s$dt > 0) "discrete" else "continuous",
                    tol = st$tol)
  ns <- sch$sdim
  T <- sch$T
  Z <- sch$Z
  Bt <- t(Z) %*% s$B
  Ct <- s$C %*% Z
  if (ns == 0L || ns == n) {
    trans <- ctrl_ss(T, Bt, Ct, s$D, dt = s$dt, input_names = s$input_names,
                     output_names = s$output_names)
    zero <- ctrl_gain(0 * s$D, dt = s$dt, input_names = s$input_names,
                      output_names = s$output_names)
    if (ns == n) return(list(stable = trans, unstable = zero))
    return(list(stable = ctrl_gain(s$D, dt = s$dt, input_names = s$input_names,
                                   output_names = s$output_names),
                unstable = ctrl_ss(T, Bt, Ct, 0 * s$D, dt = s$dt,
                                   input_names = s$input_names,
                                   output_names = s$output_names)))
  }
  i1 <- seq_len(ns)
  i2 <- ns + seq_len(n - ns)
  # T11 X - X T22 = -T12 gives [I X; 0 I]^{-1} T [I X; 0 I] block diagonal
  X <- ctrl_sylvester(T[i1, i1, drop = FALSE], -T[i2, i2, drop = FALSE],
                      -T[i1, i2, drop = FALSE])
  B1 <- Bt[i1, , drop = FALSE] - X %*% Bt[i2, , drop = FALSE]
  B2 <- Bt[i2, , drop = FALSE]
  C1 <- Ct[, i1, drop = FALSE]
  C2 <- Ct[, i1, drop = FALSE] %*% X + Ct[, i2, drop = FALSE]
  list(
    stable = ctrl_ss(T[i1, i1, drop = FALSE], B1, C1, s$D, dt = s$dt,
                     input_names = s$input_names, output_names = s$output_names),
    unstable = ctrl_ss(T[i2, i2, drop = FALSE], B2, C2, 0 * s$D, dt = s$dt,
                       input_names = s$input_names,
                       output_names = s$output_names)
  )
}

#' Hankel singular values
#'
#' The square roots of the eigenvalues of \eqn{W_c W_o} for an
#' asymptotically stable model, computed from Gramian factors as the
#' singular values of \eqn{L_o^\top L_c}.
#'
#' @param sys An asymptotically stable model.
#' @return A decreasing numeric vector.
#' @examples
#' ctrl_hsvd(ctrl_tf(1, c(1, 3, 2)))
#' @export
ctrl_hsvd <- function(sys) {
  s <- prescale_ss(as_ss(sys))
  if (!nrow(s$A)) return(numeric())
  f <- gramian_factors(s)
  if (!ncol(f$Lo) || !ncol(f$Lc)) return(rep(0, nrow(s$A)))
  hs <- svd(t(f$Lo) %*% f$Lc, nu = 0, nv = 0)$d
  c(hs, rep(0, nrow(s$A) - length(hs)))
}

gramian_factors <- function(s) {
  Wc <- ctrl_gram(s, "c")
  Wo <- ctrl_gram(s, "o")
  list(Wc = Wc, Wo = Wo, Lc = psd_sqrt_factor(Wc), Lo = psd_sqrt_factor(Wo))
}

#' Balanced realization and balanced truncation
#'
#' `ctrl_balreal()` returns a balanced realization of an asymptotically
#' stable model, in which both Gramians equal
#' \eqn{\mathrm{diag}(\sigma_1, \ldots, \sigma_r)} (the Hankel singular
#' values); states with zero Hankel singular value (to tolerance) are
#' removed. `ctrl_balred()` truncates the balanced realization to `order`
#' states using the square-root method of Tombs and Laub, which avoids
#' forming the balancing transformation from ill-conditioned Gramians.
#'
#' For an asymptotically stable model, balanced truncation preserves
#' stability and satisfies the a priori bound
#' \deqn{\|G - G_r\|_\infty \le 2 \sum_{i > r} \sigma_i}
#' (continuous and discrete time). With `method = "matchdc"` the discarded
#' states are residualized (singular perturbation approximation), which
#' matches the DC gain and satisfies the same bound.
#'
#' Unstable models are first split with [ctrl_stabsep()]; only the stable
#' part is reduced and the unstable part is kept exactly, so `order` must be
#' at least the number of unstable poles.
#'
#' @param sys A model.
#' @param order Number of states to keep.
#' @param method `"truncate"` or `"matchdc"`.
#' @param tol Relative tolerance below which Hankel singular values are
#'   treated as zero in `ctrl_balreal()`.
#' @return `ctrl_balreal()` returns a `ctrl_ss` model with attribute `hsv`.
#'   `ctrl_balred()` returns the reduced `ctrl_ss` model with attributes
#'   `hsv` (Hankel singular values of the stable part) and `error_bound`.
#' @examples
#' G <- ctrl_tf(c(1, 11, 30), c(1, 12, 47, 72, 36))
#' Gr <- ctrl_balred(G, 2)
#' attr(Gr, "error_bound")
#' ctrl_norm(ctrl_tf2ss(G) - Gr, "inf")
#' @export
ctrl_balreal <- function(sys, tol = 1e-12) {
  s <- prescale_ss(as_ss(sys))
  if (!nrow(s$A)) return(s)
  f <- gramian_factors(s)
  if (!ncol(f$Lo) || !ncol(f$Lc)) {
    # no state is both controllable and observable: the model is static
    out <- ctrl_gain(s$D, dt = s$dt, input_names = s$input_names,
                     output_names = s$output_names)
    attr(out, "hsv") <- rep(0, nrow(s$A))
    return(out)
  }
  sv <- svd(t(f$Lo) %*% f$Lc)
  r <- sum(sv$d > tol * max(sv$d))
  out <- balanced_project(s, f, sv, r, "truncate")
  attr(out, "hsv") <- sv$d
  out
}

balanced_project <- function(s, f, sv, r, method) {
  idx <- seq_len(r)
  sig <- sv$d[idx]
  Sinv <- diag(1 / sqrt(sig), r)
  Tl <- Sinv %*% t(sv$u[, idx, drop = FALSE]) %*% t(f$Lo)
  Tr <- f$Lc %*% sv$v[, idx, drop = FALSE] %*% Sinv
  rr <- sum(sv$d > 1e-14 * max(sv$d))
  if (method == "truncate" || r >= rr) {
    return(ctrl_ss(Tl %*% s$A %*% Tr, Tl %*% s$B, s$C %*% Tr, s$D, dt = s$dt,
                   input_names = s$input_names, output_names = s$output_names))
  }
  # full balancing (all states with positive Hankel singular values) then
  # residualize the trailing block
  full <- balanced_project(s, f, sv, rr, "truncate")
  i1 <- seq_len(r)
  i2 <- (r + 1L):rr
  A <- full$A
  B <- full$B
  C <- full$C
  if (s$dt > 0) {
    M <- diag(length(i2)) - A[i2, i2, drop = FALSE]
  } else {
    M <- -A[i2, i2, drop = FALSE]
  }
  X <- safe_solve(M, cbind(A[i2, i1, drop = FALSE], B[i2, , drop = FALSE]),
                  "residualized block")
  XA <- X[, seq_along(i1), drop = FALSE]
  XB <- X[, length(i1) + seq_len(ncol(B)), drop = FALSE]
  ctrl_ss(A[i1, i1, drop = FALSE] + A[i1, i2, drop = FALSE] %*% XA,
          B[i1, , drop = FALSE] + A[i1, i2, drop = FALSE] %*% XB,
          C[, i1, drop = FALSE] + C[, i2, drop = FALSE] %*% XA,
          s$D + C[, i2, drop = FALSE] %*% XB, dt = s$dt,
          input_names = s$input_names, output_names = s$output_names)
}

#' @rdname ctrl_balreal
#' @export
ctrl_balred <- function(sys, order, method = c("truncate", "matchdc")) {
  method <- match.arg(method)
  s <- as_ss(sys)
  n <- nrow(s$A)
  if (!is.numeric(order) || length(order) != 1L || order < 0 || order > n ||
      order != round(order)) {
    ctrl_abort("`order` must be an integer between 0 and ", n, ".")
  }
  if (ctrl_is_stable(s)) {
    stable <- s
    unstable <- NULL
  } else {
    parts <- ctrl_stabsep(s)
    stable <- parts$stable
    unstable <- parts$unstable
  }
  nu <- if (is.null(unstable)) 0L else nrow(unstable$A)
  if (order < nu) {
    ctrl_abort("`order` must be at least the number of unstable poles (", nu,
               ").")
  }
  rs <- order - nu
  if (nrow(stable$A)) {
    # equilibrate the state coordinates so that the Gramian factors keep
    # every direction that matters (their eigenvalues can otherwise span
    # more than the double precision range)
    stable <- prescale_ss(stable)
    f <- gramian_factors(stable)
    sv <- if (ncol(f$Lo) && ncol(f$Lc)) svd(t(f$Lo) %*% f$Lc) else list(d = numeric())
    rs <- min(rs, sum(sv$d > 1e-14 * max(sv$d, 0)))
    red <- if (rs == 0L) {
      dc <- if (method == "matchdc") ctrl_dcgain(stable) else stable$D
      ctrl_gain(dc, dt = s$dt, input_names = s$input_names,
                output_names = s$output_names)
    } else balanced_project(stable, f, sv, rs, method)
    hsv <- sv$d
  } else {
    red <- stable
    hsv <- numeric()
  }
  out <- if (is.null(unstable)) red else ctrl_parallel(red, unstable)
  out$input_names <- s$input_names
  out$output_names <- s$output_names
  bound <- 2 * sum(hsv[seq_along(hsv) > rs])
  if (nrow(stable$A)) check_reduction_bound(stable, red, bound, method)
  attr(out, "hsv") <- hsv
  attr(out, "error_bound") <- bound
  out
}

# The a priori bound must dominate the actual error. Frequency samples give
# rigorous lower bounds on the H-infinity error; if one exceeds the bound,
# the Gramian computation lost directions and the bound cannot be trusted.
check_reduction_bound <- function(full, red, bound, method) {
  w <- c(0, ctrl_freq_grid(full, 200L))
  if (full$dt > 0) w <- w[w <= pi / full$dt]
  z <- if (full$dt > 0) exp(1i * w * full$dt) else 1i * w
  gf <- eval_system(full, z)
  ef <- gf - eval_system(as_ss(red), z)
  smax <- function(v) if (any(!is.finite(v))) 0 else max(svd(as.matrix(v), 0, 0)$d, 0)
  lower <- max(vapply(seq_along(w), function(k) smax(ef[, , k]), 0))
  # evaluation error is relative to the size of the transfer matrix itself,
  # not to an absolute unit
  gscale <- max(vapply(seq_along(w), function(k) smax(gf[, , k]), 0))
  if (lower > bound + 1e-8 * gscale + 1e-6 * bound) {
    ctrl_abort("Balanced reduction could not certify its error bound: the ",
               "sampled error ", signif(lower, 4), " exceeds the bound ",
               signif(bound, 4), ". The Gramians are too ill-conditioned for ",
               "the reduction to be reliable.")
  }
  invisible(lower)
}
