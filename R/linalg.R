# Internal wrappers around LAPACK through src/linalg.c.

# Real Schur decomposition A = Z T Z'.
# `order`: "none", "continuous" (Re < -tol first) or "discrete" (|.| < 1 - tol first).
schur_real <- function(A, order = c("none", "continuous", "discrete"), tol = 0) {
  order <- match.arg(order)
  mode <- match(order, c("none", "continuous", "discrete")) - 1L
  A <- as.matrix(A)
  storage.mode(A) <- "double"
  if (!length(A)) {
    return(list(T = A, Z = A, values = complex(), sdim = 0L))
  }
  out <- .Call(cc_schur, A, mode, as.double(tol))
  if (out$info > 0L && out$info <= nrow(A)) {
    ctrl_abort("Real Schur decomposition did not converge (LAPACK dgees info ",
               out$info, ").")
  }
  if (out$info == nrow(A) + 1L || out$info == nrow(A) + 2L) {
    ctrl_abort("Reordering of the real Schur form failed (LAPACK dgees info ",
               out$info, "); eigenvalues are too close to the stability boundary.")
  }
  list(T = out$T, Z = out$Z, values = complex(real = out$wr, imaginary = out$wi),
       sdim = out$sdim)
}

# Generalized real Schur (QZ) decomposition Q' A Z = AA, Q' B Z = BB.
qz_real <- function(A, B, order = c("none", "continuous", "discrete"), tol = 0) {
  order <- match.arg(order)
  mode <- match(order, c("none", "continuous", "discrete")) - 1L
  storage.mode(A) <- "double"
  storage.mode(B) <- "double"
  out <- .Call(cc_qz, A, B, mode, as.double(tol))
  n <- nrow(A)
  if (out$info > 0L && out$info <= n) {
    ctrl_abort("QZ iteration did not converge (LAPACK dgges info ", out$info, ").")
  }
  if (out$info > n) {
    ctrl_abort("Reordering of the generalized Schur form failed (LAPACK dgges ",
               "info ", out$info, "); eigenvalues are too close to the ",
               "stability boundary or the pencil is ill-conditioned.")
  }
  list(AA = out$AA, BB = out$BB, Q = out$Q, Z = out$Z,
       alpha = complex(real = out$alphar, imaginary = out$alphai),
       beta = out$beta, sdim = out$sdim)
}

# Finite generalized eigenvalues of the pencil (A, B).
gen_eigenvalues <- function(A, B, tol = NULL, warn = TRUE) {
  if (!length(A)) return(complex())
  q <- qz_real(A, B)
  # beta is computed from an orthogonally equivalent copy of B, so its scale
  # is that of B; tiny beta means an infinite eigenvalue.
  sc <- max(norm(B, "1"), norm(A, "1"))
  tol <- tol %||% (1e3 * nrow(A) * .Machine$double.eps * (if (sc > 0) sc else 1))
  finite <- abs(q$beta) > tol
  # beta above pure rounding but inside the band: the eigenvalue may be a
  # large finite one or an infinite one; say so instead of dropping it
  # silently
  unsure <- !finite & abs(q$beta) > 10 * nrow(A) * .Machine$double.eps * sc &
    abs(q$alpha) > 0
  if (warn && any(unsure)) {
    ctrl_warn(sum(unsure), " generalized eigenvalue(s) of size about ",
              paste(signif(Mod(q$alpha[unsure] / q$beta[unsure]), 3),
                    collapse = ", "),
              " cannot be distinguished from infinity at this precision ",
              "and were treated as infinite.")
  }
  q$alpha[finite] / q$beta[finite]
}

# Solve A X + X B = C for quasi-triangular A and B (both in real Schur form).
trsyl <- function(Ta, Tb, C, trana = "N", tranb = "N", isgn = 1L) {
  storage.mode(C) <- "double"
  out <- .Call(cc_trsyl, Ta, Tb, C, trana, tranb, as.integer(isgn))
  if (out$info < 0L) ctrl_abort("Invalid argument to LAPACK dtrsyl.")
  if (out$info == 1L) {
    ctrl_warn("Sylvester equation is nearly singular; the solution was ",
              "computed with perturbed eigenvalues and may be inaccurate.")
  }
  out$X / out$scale
}

# Block structure (1 x 1 and 2 x 2 diagonal blocks) of a real Schur form.
schur_blocks <- function(T) {
  n <- nrow(T)
  blocks <- list()
  i <- 1L
  while (i <= n) {
    if (i < n && T[i + 1L, i] != 0) {
      blocks[[length(blocks) + 1L]] <- c(i, i + 1L)
      i <- i + 2L
    } else {
      blocks[[length(blocks) + 1L]] <- i
      i <- i + 1L
    }
  }
  blocks
}
