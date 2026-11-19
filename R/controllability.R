#' Controllability and observability matrices
#'
#' \deqn{\mathcal C = [B, AB, \ldots, A^{n-1}B], \qquad
#' \mathcal O = [C; CA; \ldots; CA^{n-1}].}
#' These matrices are badly scaled for all but small systems; prefer
#' [ctrl_controllability()] and [ctrl_observability()], which use an
#' orthogonal staircase reduction, for rank decisions.
#'
#' @param sys A `ctrl_ss` model, or the matrix `A` when `B` (or `C`) is given.
#' @param B,C Optional matrices when `sys` is a matrix.
#' @return A numeric matrix.
#' @examples
#' A <- matrix(c(0, 1, -2, -3), 2, byrow = TRUE)
#' ctrl_ctrb(A, c(0, 1))
#' ctrl_obsv(A, C = c(1, 0))
#' @export
ctrl_ctrb <- function(sys, B = NULL) {
  ab <- get_ab(sys, B, "B")
  A <- ab$A
  B <- ab$B
  n <- nrow(A)
  out <- matrix(0, n, n * ncol(B))
  blk <- B
  for (k in seq_len(n)) {
    out[, (k - 1L) * ncol(B) + seq_len(ncol(B))] <- blk
    blk <- A %*% blk
  }
  out
}

#' @rdname ctrl_ctrb
#' @export
ctrl_obsv <- function(sys, C = NULL) {
  if (is_ctrl_system(sys)) {
    s <- as_ss(sys)
    return(t(ctrl_ctrb(t(s$A), t(s$C))))
  }
  A <- as_real_matrix(sys, "A")
  check_square(A, "A")
  C <- as_real_matrix(C, "C", ncol = nrow(A))
  check_dims(C, "C", nrow(C), nrow(A))
  t(ctrl_ctrb(t(A), t(C)))
}

get_ab <- function(sys, B, bname) {
  if (is_ctrl_system(sys)) {
    s <- as_ss(sys)
    return(list(A = s$A, B = s$B))
  }
  A <- as_real_matrix(sys, "A")
  check_square(A, "A")
  B <- as_real_matrix(B, bname, nrow = nrow(A))
  check_dims(B, bname, nrow(A), ncol(B))
  list(A = A, B = B)
}

# Orthogonal controllability staircase (Varga). Returns U with
# U' A U = [Ac A12; 0 Anc], U' B = [Bc; 0] and the controllable dimension.
staircase <- function(A, B, tol = NULL) {
  n <- nrow(A)
  U <- diag(n)
  # rounding grows over the staircase stages, so scale with n^2. The first
  # stage is a rank decision on B and the later ones on blocks of A, so
  # each is relative to that matrix (the rank of (A, B) does not change
  # when A or B is multiplied by a constant)
  fixed <- !is.null(tol)
  rel <- 10 * max(n, ncol(B), 1)^2 * .Machine$double.eps
  tol_b <- if (fixed) tol else rel * norm(B, "2")
  tol_a <- if (fixed) tol else rel * norm(A, "2")
  tol <- tol_b
  k <- 0L
  M <- B
  while (k < n) {
    if (!length(M) || ncol(M) == 0L) break
    s <- svd(M, nu = nrow(M), nv = 0)
    r <- sum(s$d > (if (k == 0L) tol_b else tol_a))
    if (r == 0L) break
    T <- blkdiag(diag(k), s$u)
    A <- t(T) %*% A %*% T
    B <- t(T) %*% B
    U <- U %*% T
    prev <- (k + 1L):(k + r)
    k <- k + r
    if (k >= n) break
    M <- A[(k + 1L):n, prev, drop = FALSE]
  }
  list(U = U, A = A, B = B, nc = k, tol = tol)
}

pbh_table <- function(A, B, tol) {
  n <- nrow(A)
  lam <- as.complex(eigen(A, only.values = TRUE)$values)
  smin <- vapply(seq_along(lam), function(i) {
    M <- cbind(A - lam[i] * diag(n), B)
    min(svd(M, nu = 0, nv = 0)$d[seq_len(n)])
  }, 0)
  data.frame(eigenvalue = lam, pbh_sigma_min = smin,
             controllable = smin > tol)
}

#' Controllability and stabilizability analysis
#'
#' Determine the controllable subspace with an orthogonal staircase
#' reduction, which is far better conditioned than the rank of the
#' controllability matrix. The uncontrollable modes are the eigenvalues of
#' the uncontrollable block. The model is *stabilizable* when every
#' uncontrollable mode is asymptotically stable; full controllability is not
#' required.
#'
#' Rank decisions use the tolerance `tol`, scaled by the norms of `A` and
#' `B`. A rank computed under a tolerance is a numerical statement: a
#' nearly uncontrollable system can be reported as controllable with a
#' tiny PBH singular value. Inspect `pbh$pbh_sigma_min` (the smallest
#' singular value of \eqn{[A - \lambda I, B]} for each eigenvalue) to judge
#' the distance to uncontrollability.
#'
#' @param sys A `ctrl_ss` model (or the matrix `A`).
#' @param B Input matrix when `sys` is a matrix.
#' @param tol Rank tolerance.
#' @param discrete Time domain when `sys` is a matrix (`FALSE` by default).
#' @return A list of class `ctrl_controllability` with `rank`, `n`,
#'   `controllable`, `stabilizable`, `uncontrollable_modes`, `basis` (an
#'   orthonormal basis of the controllable subspace), `transform` (the
#'   orthogonal staircase transformation), `pbh` (table of PBH tests) and
#'   `ctrb_singular_values`.
#' @examples
#' A <- diag(c(-1, 2))
#' B <- c(1, 0)
#' res <- ctrl_controllability(A, B)
#' res$controllable
#' res$stabilizable  # the unstable mode at 2 is uncontrollable
#' @export
ctrl_controllability <- function(sys, B = NULL, tol = NULL, discrete = NULL) {
  ab <- get_ab(sys, B, "B")
  dt <- if (is_ctrl_system(sys)) sys$dt else if (isTRUE(discrete)) 1 else 0
  controllability_core(ab$A, ab$B, tol, dt, "controllable")
}

controllability_core <- function(A, B, tol, dt, word) {
  n <- nrow(A)
  st <- staircase(A, B, tol)
  nc <- st$nc
  unc <- if (nc < n) {
    Anc <- st$A[(nc + 1L):n, (nc + 1L):n, drop = FALSE]
    as.complex(eigen(Anc, only.values = TRUE)$values)
  } else complex()
  stab_tol <- default_eig_tol(A)
  stable_mode <- if (dt > 0) Mod(unc) < 1 - stab_tol else Re(unc) < -stab_tol
  out <- list(
    rank = nc,
    n = n,
    full = nc == n,
    stabilizable = all(stable_mode),
    uncontrollable_modes = unc,
    basis = st$U[, seq_len(nc), drop = FALSE],
    transform = st$U,
    pbh = if (n) pbh_table(A, B, st$tol) else NULL,
    ctrb_singular_values = if (n) ctrb_sv(A, B) else numeric(),
    tol = st$tol,
    kind = word
  )
  if (word == "controllable") {
    out$controllable <- out$full
  } else {
    names(out)[names(out) == "stabilizable"] <- "detectable"
    names(out)[names(out) == "uncontrollable_modes"] <- "unobservable_modes"
    out$observable <- out$full
    out$pbh$observable <- out$pbh$controllable
    out$pbh$controllable <- NULL
  }
  class(out) <- "ctrl_controllability"
  out
}

#' Observability and detectability analysis
#'
#' The dual of [ctrl_controllability()]: the analysis is applied to
#' \eqn{(A^\top, C^\top)}. The model is *detectable* when every unobservable
#' mode is asymptotically stable.
#'
#' @inheritParams ctrl_controllability
#' @param C Output matrix when `sys` is a matrix.
#' @return A list of class `ctrl_controllability` with `rank`, `observable`,
#'   `detectable`, `unobservable_modes`, `basis` (orthonormal basis of the
#'   observable subspace, that is the orthogonal complement of the
#'   unobservable subspace), `pbh` and `ctrb_singular_values` (singular
#'   values of the observability matrix).
#' @examples
#' ctrl_observability(diag(c(-1, -2)), C = c(1, 0))$detectable
#' @export
ctrl_observability <- function(sys, C = NULL, tol = NULL, discrete = NULL) {
  if (is_ctrl_system(sys)) {
    s <- as_ss(sys)
    A <- s$A
    C <- s$C
    dt <- s$dt
  } else {
    A <- as_real_matrix(sys, "A")
    check_square(A, "A")
    C <- as_real_matrix(C, "C", ncol = nrow(A))
    dt <- if (isTRUE(discrete)) 1 else 0
  }
  controllability_core(t(A), t(C), tol, dt, "observable")
}

#' @export
print.ctrl_controllability <- function(x, digits = 4, ...) {
  obs <- x$kind == "observable"
  cat(if (obs) "Observability" else "Controllability", " analysis: rank ",
      x$rank, " of ", x$n, " (tol = ", signif(x$tol, 3), ")\n", sep = "")
  cat("  ", if (obs) "observable" else "controllable", ": ", x$full, "\n",
      sep = "")
  cat("  ", if (obs) "detectable" else "stabilizable", ": ",
      if (obs) x$detectable else x$stabilizable, "\n", sep = "")
  modes <- if (obs) x$unobservable_modes else x$uncontrollable_modes
  if (length(modes)) {
    cat("  ", if (obs) "unobservable" else "uncontrollable", " modes: ",
        paste(format(signif(modes, digits)), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

# Singular values of the controllability matrix for diagnostics only; NA
# when the matrix powers overflow (the staircase result does not need them).
ctrb_sv <- function(A, B) {
  M <- ctrl_ctrb(A, B)
  if (any(!is.finite(M))) return(rep(NA_real_, nrow(A)))
  tryCatch(svd(M, 0, 0)$d, error = function(e) rep(NA_real_, nrow(A)))
}

#' Stabilizability and detectability tests
#'
#' @inheritParams ctrl_controllability
#' @return `TRUE` or `FALSE`.
#' @examples
#' ctrl_is_stabilizable(ctrl_ss(diag(c(-1, 2)), c(0, 1), diag(2)))
#' ctrl_is_detectable(ctrl_ss(diag(c(-1, 2)), c(0, 1), c(1, 0)))
#' @export
ctrl_is_stabilizable <- function(sys, tol = NULL) {
  ctrl_controllability(sys, tol = tol)$stabilizable
}

#' @rdname ctrl_is_stabilizable
#' @export
ctrl_is_detectable <- function(sys, tol = NULL) {
  ctrl_observability(sys, tol = tol)$detectable
}

#' @rdname ctrl_is_stabilizable
#' @export
ctrl_is_controllable <- function(sys, tol = NULL) {
  ctrl_controllability(sys, tol = tol)$full
}

#' @rdname ctrl_is_stabilizable
#' @export
ctrl_is_observable <- function(sys, tol = NULL) {
  ctrl_observability(sys, tol = tol)$full
}
