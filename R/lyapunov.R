#' Lyapunov and Stein equations
#'
#' Solve the continuous time Lyapunov equation
#' \deqn{A X + X A^\top + Q = 0}
#' or the discrete time Lyapunov (Stein) equation
#' \deqn{A X A^\top - X + Q = 0.}
#'
#' Both solvers reduce `A` to real Schur form with LAPACK. The continuous
#' equation is then solved with the Bartels-Stewart algorithm (LAPACK
#' `dtrsyl`); the discrete equation with a block back substitution over the
#' quasi-triangular Schur factor (Barraud's method). No Kronecker
#' vectorization of the full problem is used, so the cost is
#' \eqn{O(n^3)}.
#'
#' The solution is unique exactly when \eqn{\lambda_i + \lambda_j \ne 0}
#' (continuous) or \eqn{\lambda_i \lambda_j \ne 1} (discrete) for all pairs
#' of eigenvalues of `A`. If the smallest such separation, relative to the
#' scale of `A`, is below `tol`, an error is raised instead of returning an
#' unreliable solution.
#'
#' @param A Square matrix.
#' @param Q Square matrix of the same size (usually symmetric).
#' @param discrete `FALSE` (default) for the continuous equation, `TRUE` for
#'   the discrete equation.
#' @param tol Relative uniqueness tolerance.
#' @return The solution matrix `X`, with attributes `residual` (relative
#'   Frobenius residual of the equation) and `separation` (the uniqueness
#'   measure). When `Q` is symmetric the solution is symmetrized.
#' @examples
#' A <- matrix(c(-1, 2, 0, -3), 2)
#' X <- ctrl_lyap(A, diag(2))
#' A %*% X + X %*% t(A) + diag(2)  # close to zero
#' Xd <- ctrl_lyap(A / 4, diag(2), discrete = TRUE)
#' attr(Xd, "residual")
#' @export
ctrl_lyap <- function(A, Q, discrete = FALSE, tol = 1e-12) {
  A <- as_real_matrix(A, "A")
  check_square(A, "A")
  n <- nrow(A)
  Q <- as_real_matrix(Q, "Q", n, n)
  check_dims(Q, "Q", n, n)
  check_flag(discrete, "discrete")
  if (n == 0L) return(matrix(0, 0, 0))
  sym_q <- asymmetry(Q) <= 1e-12
  sch <- schur_real(A)
  lam <- sch$values
  scale <- max(1, max(Mod(lam)))
  sep <- if (discrete) {
    min(Mod(1 - outer(lam, lam))) / scale^2
  } else {
    min(Mod(outer(lam, lam, "+"))) / scale
  }
  if (sep < tol) {
    ctrl_abort("The ", if (discrete) "discrete" else "continuous",
               " Lyapunov equation has no unique solution: A has eigenvalues ",
               if (discrete) "with lambda_i * lambda_j = 1"
               else "with lambda_i + lambda_j = 0",
               " (relative separation ", signif(sep, 3), ").")
  }
  Z <- sch$Z
  T <- sch$T
  Qt <- t(Z) %*% Q %*% Z
  Xt <- if (discrete) stein_schur(T, Qt) else trsyl(T, T, -Qt, "N", "T", 1L)
  X <- Z %*% Xt %*% t(Z)
  if (sym_q) X <- symmetrize(X)
  R <- if (discrete) A %*% X %*% t(A) - X + Q else A %*% X + X %*% t(A) + Q
  attr(X, "residual") <- rel_residual(R, Q, A %*% X)
  attr(X, "separation") <- sep
  X
}

# Solve T X T' - X + Q = 0 for quasi upper triangular T (real Schur form).
stein_schur <- function(T, Q) {
  n <- nrow(T)
  X <- matrix(0, n, n)
  blocks <- schur_blocks(T)
  p <- length(blocks)
  for (jb in p:1) {
    J <- blocks[[jb]]
    later <- if (max(J) < n) (max(J) + 1L):n else integer()
    # contribution of already solved columns: X[, later] T[J, later]'
    Ypre <- if (length(later)) {
      X[, later, drop = FALSE] %*% t(T[J, later, drop = FALSE])
    } else matrix(0, n, length(J))
    Tjj <- T[J, J, drop = FALSE]
    for (ib in p:1) {
      I <- blocks[[ib]]
      rows_ge <- min(I):n
      rows_gt <- if (max(I) < n) (max(I) + 1L):n else integer()
      rhs <- -Q[I, J, drop = FALSE] -
        T[I, rows_ge, drop = FALSE] %*% Ypre[rows_ge, , drop = FALSE]
      if (length(rows_gt)) {
        rhs <- rhs - T[I, rows_gt, drop = FALSE] %*%
          X[rows_gt, J, drop = FALSE] %*% t(Tjj)
      }
      Tii <- T[I, I, drop = FALSE]
      # Tii Xij Tjj' - Xij = rhs  <=>  (Tjj kron Tii - I) vec(Xij) = vec(rhs)
      K <- kronecker(Tjj, Tii) - diag(length(I) * length(J))
      X[I, J] <- matrix(solve(K, as.vector(rhs)), length(I), length(J))
    }
  }
  X
}

#' Sylvester equation
#'
#' Solve \eqn{A X + X B = C} with the Bartels-Stewart algorithm. A unique
#' solution exists when `A` and `-B` have no common eigenvalues.
#'
#' @param A Square `n x n` matrix.
#' @param B Square `m x m` matrix.
#' @param C `n x m` matrix.
#' @param tol Relative uniqueness tolerance.
#' @return The `n x m` solution with a `residual` attribute.
#' @examples
#' X <- ctrl_sylvester(diag(c(1, 2)), matrix(3), matrix(1, 2, 1))
#' X
#' @export
ctrl_sylvester <- function(A, B, C, tol = 1e-12) {
  A <- as_real_matrix(A, "A")
  B <- as_real_matrix(B, "B")
  check_square(A, "A")
  check_square(B, "B")
  C <- as_real_matrix(C, "C", nrow(A), nrow(B))
  check_dims(C, "C", nrow(A), nrow(B))
  if (!length(C)) return(C)
  sa <- schur_real(A)
  sb <- schur_real(B)
  scale <- max(1, Mod(sa$values), Mod(sb$values))
  sep <- min(Mod(outer(sa$values, sb$values, "+"))) / scale
  if (sep < tol) {
    ctrl_abort("The Sylvester equation has no unique solution: A and -B ",
               "share an eigenvalue (relative separation ", signif(sep, 3), ").")
  }
  Ct <- t(sa$Z) %*% C %*% sb$Z
  Xt <- trsyl(sa$T, sb$T, Ct, "N", "N", 1L)
  X <- sa$Z %*% Xt %*% t(sb$Z)
  attr(X, "residual") <- rel_residual(A %*% X + X %*% B - C, C, A %*% X)
  X
}

#' Controllability and observability Gramians
#'
#' For an asymptotically stable model, the controllability Gramian solves
#' \eqn{A W_c + W_c A^\top + B B^\top = 0} (continuous) or
#' \eqn{A W_c A^\top - W_c + B B^\top = 0} (discrete), and the
#' observability Gramian solves the dual equations with
#' \eqn{A^\top} and \eqn{C^\top C}.
#'
#' @param sys A `ctrl_ss` model (asymptotically stable).
#' @param type `"c"` (controllability) or `"o"` (observability).
#' @return A symmetric positive semidefinite matrix.
#' @examples
#' sys <- ctrl_ss(matrix(c(-1, 0, 1, -2), 2), c(0, 1), c(1, 0))
#' ctrl_gram(sys, "c")
#' @export
ctrl_gram <- function(sys, type = c("c", "o")) {
  type <- match.arg(type)
  sys <- as_ss(sys)
  if (!ctrl_is_stable(sys)) {
    ctrl_abort("Gramians are only defined here for asymptotically stable ",
               "models; separate the unstable part first with ctrl_stabsep().")
  }
  d <- sys$dt > 0
  if (type == "c") {
    ctrl_lyap(sys$A, sys$B %*% t(sys$B), discrete = d)
  } else {
    ctrl_lyap(t(sys$A), t(sys$C) %*% sys$C, discrete = d)
  }
}
