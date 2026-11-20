#' Poles of a model
#'
#' For a state-space model the poles are the eigenvalues of `A`. For a
#' transfer function they are the roots of the denominator polynomials
#' (after the coefficient normalization done by [ctrl_tf()]); pole-zero
#' cancellations are not removed.
#'
#' @param sys A `ctrl_ss`, `ctrl_tf` or `ctrl_zpk` model.
#' @return A complex vector of poles sorted by decreasing real part
#'   (continuous time) or decreasing modulus (discrete time).
#' @examples
#' ctrl_poles(ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#' ctrl_poles(ctrl_tf(1, c(1, 3, 2)))
#' @export
ctrl_poles <- function(sys) {
  assert_system(sys)
  p <- if (is_ctrl_tf(sys)) {
    unlist(lapply(unique_denominators(sys), poly_roots))
  } else if (inherits(sys, "ctrl_zpk")) {
    unlist(sys$poles)
  } else {
    if (!nrow(sys$A)) return(complex())
    as.complex(eigen(sys$A, only.values = TRUE)$values)
  }
  sort_poles(as.complex(p), sys$dt)
}

sort_poles <- function(p, dt) {
  if (!length(p)) return(complex())
  key <- if (dt > 0) Mod(p) else Re(p)
  p[order(-key, -Im(p))]
}

#' Stability analysis
#'
#' Classify the internal stability of a linear time-invariant model from its
#' poles, with an explicit numerical tolerance.
#'
#' A continuous time model is asymptotically stable when every pole has
#' real part below `-tol`; a discrete time model when every pole has modulus
#' below `1 - tol`. Poles within `tol` of the stability boundary (the
#' imaginary axis or the unit circle) are *boundary poles*. If there are
#' boundary poles but no poles outside the boundary, the model is
#' *marginally stable* when every boundary eigenvalue is semisimple
#' (geometric multiplicity equals algebraic multiplicity, so trajectories
#' remain bounded) and *unstable* otherwise (a Jordan block on the boundary
#' produces polynomial growth). Multiplicities are numerical estimates:
#' eigenvalues are clustered within `cluster_tol` and the geometric
#' multiplicity uses a scale-aware SVD rank.
#'
#' For transfer functions only the denominator roots are available, so the
#' Jordan structure test cannot be done; repeated boundary poles are then
#' reported as unstable (true for the input-output response).
#'
#' @param sys A model.
#' @param tol Boundary tolerance. Defaults to a multiple of machine precision
#'   times the norm of the balanced `A` (the state coordinates are first
#'   scaled to balance \eqn{[A\ B; C\ 0]}, so the result does not depend on
#'   state units, and a uniformly slow continuous model keeps its relative
#'   accuracy; in discrete time the scale is at least one). For a model
#'   whose time scales differ by about 1e13 or more, the slowest stable
#'   poles can fall inside this band and are classified as marginal; pass a
#'   smaller `tol` when such poles are known to be accurate.
#' @param cluster_tol Tolerance for grouping eigenvalues into clusters.
#' @return An object of class `ctrl_stability` with `classification`
#'   (`"asymptotically stable"`, `"marginally stable"` or `"unstable"`),
#'   `stable` (asymptotic stability), `poles`, `boundary_poles`,
#'   `unstable_poles`, `stability_margin` (distance of the slowest pole to
#'   the boundary; negative when unstable), the uniform band `tol`, the
#'   per-pole bands `pole_tol` actually used and `cluster_tol`.
#' @examples
#' ctrl_stability(ctrl_ss(-1, 1, 1))
#' # A double integrator is not marginally stable (Jordan block at 0)
#' ctrl_stability(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#' # A harmonic oscillator is marginally stable
#' ctrl_stability(ctrl_ss(matrix(c(0, 1, -1, 0), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#' @export
ctrl_stability <- function(sys, tol = NULL, cluster_tol = NULL) {
  assert_system(sys)
  discrete <- sys$dt > 0
  A <- NULL
  if (is_ctrl_ss(sys)) {
    # a diagonal similarity that balances [A B; C 0] removes state-unit
    # effects (also for triangular A, which LAPACK balancing of A alone
    # cannot scale); the eigenvalues are computed from the balanced A
    A <- prescale_ss(sys)$A
    if (nrow(A)) {
      A <- A * outer(1 / (s <- expm::balance(A, job = "S")$scale), s)
    }
    poles <- complex()
    kappa <- numeric()
    if (nrow(A)) {
      ev <- eigen(A)
      V <- ev$vectors
      W <- tryCatch(solve(V), error = function(e) NULL)
      # eigenvalue condition numbers ||w_i|| ||v_i|| / |w_i v_i| (W V = I);
      # a defective (or nearly defective) matrix gets Inf
      kappa <- if (is.null(W)) rep(Inf, nrow(A)) else
        sqrt(colSums(Mod(V)^2)) * sqrt(rowSums(Mod(W)^2))
      p <- as.complex(ev$values)
      key <- if (sys$dt > 0) Mod(p) else Re(p)
      ord <- order(-key, -Im(p))
      poles <- p[ord]
      kappa <- kappa[ord]
    }
  } else {
    poles <- ctrl_poles(sys)
    kappa <- rep(Inf, length(poles))
  }
  # rounding errors in the poles scale with the balanced matrix norm (or
  # the pole size), with no fixed unit floor in continuous time: a uniformly
  # slow model keeps its relative accuracy. The unit circle sets an absolute
  # scale in discrete time.
  scale <- if (!is.null(A)) norm(A, "1") else max(0, Mod(poles))
  if (discrete) scale <- max(1, scale)
  user_tol <- !is.null(tol)
  tol <- tol %||% (1e3 * max(1, length(poles)) * .Machine$double.eps * scale)
  cluster_tol <- cluster_tol %||% max(sqrt(.Machine$double.eps) * scale, tol)
  # per-pole band: the first-order error bound of a backward stable
  # eigensolver, eps ||A|| kappa_i, when it is tighter than the uniform
  # band (well-conditioned poles of a stiff model are then classified)
  pole_tol <- if (user_tol || is.null(A)) rep(tol, length(poles)) else
    pmin(tol, 10 * max(1, length(poles)) * .Machine$double.eps *
           norm(A, "1") * kappa)
  dist <- if (discrete) 1 - Mod(poles) else -Re(poles)
  boundary <- abs(dist) <= pole_tol
  outside <- dist < -pole_tol
  classification <- "asymptotically stable"
  defective <- FALSE
  if (any(outside)) {
    classification <- "unstable"
  } else if (any(boundary)) {
    defective <- boundary_defective(poles[boundary], poles, A, cluster_tol)
    classification <- if (defective) "unstable" else "marginally stable"
  }
  out <- list(
    classification = classification,
    stable = classification == "asymptotically stable",
    poles = poles,
    boundary_poles = poles[boundary],
    unstable_poles = poles[outside],
    defective_boundary = defective,
    stability_margin = if (length(dist)) min(dist) else Inf,
    discrete = discrete,
    tol = tol,
    pole_tol = pole_tol,
    cluster_tol = cluster_tol
  )
  class(out) <- "ctrl_stability"
  out
}

boundary_defective <- function(bpoles, poles, A, cluster_tol) {
  remaining <- bpoles
  while (length(remaining)) {
    lam <- remaining[1]
    in_cluster <- Mod(poles - lam) <= cluster_tol
    alg <- sum(in_cluster)
    remaining <- remaining[Mod(remaining - lam) > cluster_tol]
    if (alg <= 1L) next
    if (is.null(A)) return(TRUE)
    n <- nrow(A)
    M <- A - lam * diag(n)
    sv <- svd(M, nu = 0, nv = 0)$d
    geo <- sum(sv <= max(cluster_tol, n * .Machine$double.eps * max(sv)))
    if (geo < alg) return(TRUE)
  }
  FALSE
}

#' @export
print.ctrl_stability <- function(x, digits = 4, ...) {
  cat("Stability: ", x$classification, " (",
      if (x$discrete) "discrete" else "continuous", " time)\n", sep = "")
  if (length(x$poles)) {
    cat("Poles: ", paste(format(signif(x$poles, digits)), collapse = ", "),
        "\n", sep = "")
  }
  if (length(x$boundary_poles)) {
    cat("Boundary poles (within tol = ", signif(x$tol, 3), "): ",
        paste(format(signif(x$boundary_poles, digits)), collapse = ", "),
        "\n", sep = "")
  }
  if (length(x$unstable_poles)) {
    cat("Unstable poles: ",
        paste(format(signif(x$unstable_poles, digits)), collapse = ", "),
        "\n", sep = "")
  }
  cat("Stability margin: ", signif(x$stability_margin, digits), "\n", sep = "")
  invisible(x)
}

#' Test asymptotic stability
#'
#' @inheritParams ctrl_stability
#' @return `TRUE` when the model is asymptotically stable.
#' @examples
#' ctrl_is_stable(ctrl_ss(-1, 1, 1))
#' @export
ctrl_is_stable <- function(sys, tol = NULL) ctrl_stability(sys, tol = tol)$stable

#' Natural frequencies and damping ratios
#'
#' For a continuous time pole \eqn{\lambda}, the natural frequency is
#' \eqn{\omega_n = |\lambda|} and the damping ratio is
#' \eqn{\zeta = -\mathrm{Re}(\lambda)/|\lambda|}. Discrete time poles are
#' first mapped to the equivalent continuous time pole
#' \eqn{\log(\lambda)/dt}; poles at the origin of the z plane have no
#' continuous equivalent and get `NA` values.
#'
#' @param sys A model.
#' @return A data frame with `pole`, `damping`, `frequency` (rad per time
#'   unit) and `time_constant` (`1 / (zeta * wn)`, `Inf` when undamped).
#' @examples
#' ctrl_damp(ctrl_tf(1, c(1, 0.4, 4)))
#' @export
ctrl_damp <- function(sys) {
  p <- ctrl_poles(sys)
  s <- p
  if (sys$dt > 0) {
    s <- ifelse(Mod(p) > 0, log(p) / sys$dt, NA_complex_)
  }
  wn <- Mod(s)
  zeta <- ifelse(wn > 0, -Re(s) / wn, -1)
  zeta[wn == 0 & !is.na(wn)] <- NA_real_
  tc <- ifelse(is.na(zeta) | is.na(wn), NA_real_,
               ifelse(zeta * wn > 0, 1 / (zeta * wn), Inf))
  data.frame(pole = p, damping = zeta, frequency = wn, time_constant = tc)
}

#' Steady-state (DC) gain
#'
#' The gain at zero frequency: \eqn{D - C A^{-1} B} in continuous time and
#' \eqn{D + C (I - A)^{-1} B} in discrete time. An error is raised when the
#' input-output map has a pole at `s = 0` (or `z = 1`), because the gain is
#' then infinite or undefined; poles there that are uncontrollable or
#' unobservable are removed first.
#'
#' @param sys A model.
#' @return A `p x m` matrix.
#' @examples
#' ctrl_dcgain(ctrl_tf(2, c(1, 3, 2)))
#' @export
ctrl_dcgain <- function(sys) {
  # balancing makes the pole-at-origin test independent of state units
  s0 <- as_ss(sys)
  s <- prescale_ss(s0)
  s$output_names <- s0$output_names
  s$input_names <- s0$input_names
  n <- nrow(s$A)
  if (!n) return(s$D)
  M <- if (s$dt > 0) diag(n) - s$A else -s$A
  rc <- rcond(M)
  if (rc < n * .Machine$double.eps && n > 0L) {
    # a pole at the origin may be uncontrollable or unobservable and then
    # does not affect the gain; retry on a minimal realization
    r <- ctrl_minreal(s)
    if (nrow(r$A) < n) {
      out <- ctrl_dcgain(r)
      dimnames(out) <- list(s$output_names, s$input_names)
      return(out)
    }
    ctrl_abort("The DC gain is undefined: the model has a pole at ",
               if (s$dt > 0) "z = 1" else "s = 0", " (rcond ", signif(rc, 3), ").")
  }
  out <- s$D + s$C %*% solve(M, s$B)
  dimnames(out) <- list(s$output_names, s$input_names)
  out
}

#' Invariant zeros
#'
#' For a square state-space model the invariant zeros are the finite
#' generalized eigenvalues of the Rosenbrock system pencil
#' \deqn{\begin{pmatrix} A & B \\ C & D \end{pmatrix} - \lambda
#' \begin{pmatrix} I & 0 \\ 0 & 0 \end{pmatrix},}
#' computed with the QZ algorithm. When `D` is well-conditioned relative to
#' its own size (reciprocal condition number above \eqn{\sqrt\epsilon}),
#' there are exactly `n` finite zeros; zeros too large for the pencil to
#' separate from infinity (from a tiny feedthrough) are then taken as the
#' largest eigenvalues of \eqn{A - B D^{-1} C}, so they are kept. Otherwise a generalized
#' eigenvalue that cannot be told apart from infinity at working precision
#' is treated as infinite with a warning. Zeros need a regular pencil (full
#' normal rank); this is checked by evaluating the transfer matrix at
#' points away from every computed eigenvalue and pole, and a singular
#' pencil raises an error. A SISO transfer function returns its
#' numerator roots. MIMO transfer functions are converted to state space,
#' so they must be square as well; nonsquare models are not supported in
#' this version.
#'
#' @param sys A model.
#' @return A complex vector of zeros.
#' @examples
#' ctrl_zeros(ctrl_tf(c(1, 2), c(1, 3, 2)))
#' @export
ctrl_zeros <- function(sys) {
  assert_system(sys)
  if (is_ctrl_tf(sys)) {
    if (nrow(sys$num) == 1L && ncol(sys$num) == 1L) {
      return(poly_roots(sys$num[[1]]))
    }
    sys <- ctrl_tf2ss(sys)
  }
  if (inherits(sys, "ctrl_zpk")) sys <- as_ss(sys)
  n <- nrow(sys$A)
  m <- ncol(sys$B)
  p <- nrow(sys$C)
  if (p != m) {
    ctrl_abort("Invariant zeros are only implemented for square models ",
               "(equal numbers of inputs and outputs).")
  }
  if (!n) return(complex())
  # with a well-conditioned feedthrough (relative to its own size, whatever
  # its units) the pencil has exactly m infinite eigenvalues, so there are
  # exactly n finite zeros: those of A - B D^-1 C
  regular_D <- rcond(sys$D) > sqrt(.Machine$double.eps)
  M <- rbind(cbind(sys$A, sys$B), cbind(sys$C, sys$D))
  N <- blkdiag(diag(n), matrix(0, m, m))
  # A diagonal similarity leaves N unchanged, so balancing M equilibrates
  # the pencil without changing its generalized eigenvalues. This makes the
  # result independent of the state coordinates' units.
  M <- expm::balance(M, job = "S")$z
  # scale N to the size of M so that finite and infinite eigenvalues are
  # separated on a common scale, then restore the units
  c <- max(norm(M, "1"), .Machine$double.xmin)
  z <- gen_eigenvalues(M, c * N, warn = !regular_D) * c
  if (regular_D && length(z) < n) {
    # zeros so large that the pencil cannot tell them from infinity: take
    # them from A - B D^-1 C, whose largest eigenvalues are accurate
    # relative to their size (the QZ values stay for the moderate zeros,
    # which A - B D^-1 C would resolve only to eps ||B D^-1 C||)
    zd <- as.complex(eigen(sys$A - sys$B %*% solve(sys$D, sys$C),
                           only.values = TRUE)$values)
    z <- c(z, zd[order(Mod(zd), decreasing = TRUE)][seq_len(n - length(z))])
  }
  # Zeros are only defined when the pencil is regular, i.e. the transfer
  # matrix has full normal rank. If the pencil is regular, G(s) has full
  # rank at every point that is neither a zero nor a pole, and every zero
  # is among the computed eigenvalues z; if it is singular, G(s) is rank
  # deficient everywhere. So one full-rank evaluation at a point away from
  # all computed eigenvalues and poles proves regularity, and rank
  # deficiency at several such points shows a singular pencil. Ranks are
  # tested on G(s) itself on the scale of the poles (the rank of the
  # pencil M - s N is badly conditioned for fast or slow models).
  pe <- eigen(sys$A, only.values = TRUE)$values
  pm <- Mod(pe)
  pm <- pm[pm > 0]
  sc <- if (length(pm)) exp(mean(log(range(pm)))) else 1
  avoid <- c(z[is.finite(z)], pe)
  # candidate points spread on circles around the origin (golden angle),
  # in the s plane scaled by sc, mapped to the z plane in discrete time
  k <- 0:63
  cand <- sc * (0.5 + (k %% 8) / 4) * exp(1i * (0.3 + 2.39996323 * k))
  if (sys$dt > 0) cand <- exp(cand / sc)
  far <- vapply(cand, function(s) {
    !length(avoid) || min(Mod(avoid - s)) > 1e-2 * max(Mod(s), 1e-300)
  }, logical(1))
  s0 <- cand[far][seq_len(min(4L, sum(far)))]
  full <- vapply(s0, function(s) {
    G0 <- eval_system(sys, s)[, , 1]
    all(is.finite(G0)) && numerical_rank(as.matrix(G0))$rank == m
  }, logical(1))
  if (!any(full)) {
    ctrl_abort("The system pencil is singular (the transfer matrix does not ",
               "have full normal rank), so its invariant zeros are not ",
               "defined. This happens, for example, when an input or output ",
               "is disconnected.")
  }
  sort_poles(z, sys$dt)
}
