#' Algebraic Riccati equations
#'
#' `ctrl_care()` computes the stabilizing solution of the continuous time
#' algebraic Riccati equation
#' \deqn{A^\top P + P A - (P B + S) R^{-1} (B^\top P + S^\top) + Q = 0,}
#' and `ctrl_dare()` the stabilizing solution of the discrete time equation
#' \deqn{A^\top P A - P - (A^\top P B + S)(R + B^\top P B)^{-1}
#' (B^\top P A + S^\top) + Q = 0.}
#'
#' Both solvers use the extended (Arnold and Laub) matrix pencil of
#' dimension \eqn{2n + m}, compress away the input block with an orthogonal
#' QR factorization (so `R` is never inverted explicitly), and compute an
#' ordered generalized real Schur decomposition with LAPACK `dgges` that
#' places the `n` stable eigenvalues first. The solution is
#' \eqn{P = U_{21} U_{11}^{-1}}, obtained with a linear solve. The discrete
#' pencil formulation does not require `A` to be invertible.
#'
#' The stabilizing solution exists and is unique when `R` is positive
#' definite, \eqn{(A, B)} is stabilizable and the pencil has no eigenvalues
#' on the stability boundary (guaranteed, for example, when
#' \eqn{Q - S R^{-1} S^\top \succeq 0} and the corresponding pair is
#' detectable). When the boundary condition fails, or when \eqn{U_{11}} is
#' singular, an error is raised rather than an invalid solution returned.
#' A solution whose relative residual exceeds 1e-6, or whose closed loop is
#' not stable, raises an error. The returned object reports the relative
#' residual, the symmetry error and the closed-loop eigenvalues so that the
#' result can be verified. Before solving, the equation is equilibrated: a
#' symplectic diagonal state scaling from balancing the Hamiltonian matrix
#' (so states in very different units carry comparable weight), a common
#' time scale (continuous time) and a common weight factor, all by powers
#' of two. The residual is reported for the equilibrated equation. The
#' Schur method is backward stable in the normwise sense, so weights that
#' differ by many orders of magnitude within one state or input direction
#' can still lose relative accuracy.
#'
#' @param A State matrix (`n x n`).
#' @param B Input matrix (`n x m`).
#' @param Q State weight (`n x n`, symmetric).
#' @param R Input weight (`m x m`, symmetric positive definite).
#' @param S Optional cross weight (`n x m`); defaults to zero.
#' @param check_weights If `TRUE` (default), require
#'   \eqn{\begin{pmatrix} Q & S \\ S^\top & R\end{pmatrix} \succeq 0}, the
#'   condition under which the equation describes a well-posed LQ problem.
#'   Set to `FALSE` to solve indefinite equations (for example in
#'   H-infinity synthesis).
#' @return A list of class `ctrl_riccati` with `P` (stabilizing solution),
#'   `K` (the associated gain \eqn{R^{-1}(B^\top P + S^\top)} or
#'   \eqn{(R + B^\top P B)^{-1}(B^\top P A + S^\top)}), `poles`
#'   (eigenvalues of \eqn{A - B K}), `residual` (relative Frobenius
#'   residual of the equilibrated equation), `symmetry_error`, `rcond_u11` and `discrete`.
#' @references Arnold, W. F. and Laub, A. J. (1984). Generalized eigenproblem
#'   algorithms and software for algebraic Riccati equations. *Proceedings of
#'   the IEEE*, 72(12), 1746-1754.
#' @examples
#' # Scalar CARE with known solution sqrt(2) - 1
#' ctrl_care(-1, 1, 1, 1)$P
#' # Scalar DARE with known solution (1 + sqrt(5)) / 2
#' ctrl_dare(1, 1, 1, 1)$P
#' @export
ctrl_care <- function(A, B, Q, R, S = NULL, check_weights = TRUE) {
  riccati_solve(A, B, Q, R, S, discrete = FALSE, check_weights = check_weights)
}

#' @rdname ctrl_care
#' @export
ctrl_dare <- function(A, B, Q, R, S = NULL, check_weights = TRUE) {
  riccati_solve(A, B, Q, R, S, discrete = TRUE, check_weights = check_weights)
}

riccati_inputs <- function(A, B, Q, R, S, check_weights) {
  A <- as_real_matrix(A, "A")
  check_square(A, "A")
  n <- nrow(A)
  B <- as_real_matrix(B, "B", nrow = n)
  check_dims(B, "B", n, ncol(B))
  m <- ncol(B)
  Q <- check_symmetric(as_real_matrix(Q, "Q", n, n), "Q")
  check_dims(Q, "Q", n, n)
  R <- check_symmetric(as_real_matrix(R, "R", m, m), "R")
  check_dims(R, "R", m, m)
  S <- if (is.null(S)) matrix(0, n, m) else as_real_matrix(S, "S", n, m)
  check_dims(S, "S", n, m)
  if (check_weights) {
    check_pd(R, "R")
    # with R positive definite, [Q S; S' R] is positive semidefinite exactly
    # when the Schur complement Q - S R^-1 S' is; testing it on its own scale
    # keeps a large R from masking an indefinite Q
    SRS <- S %*% solve(R, t(S))
    Sc <- symmetrize(Q - SRS)
    tolw <- 1e3 * n * .Machine$double.eps * max(norm(Q, "2"), norm(SRS, "2"))
    if (min_eig_sym(Sc) < -tolw) {
      ctrl_abort("The joint weight [Q S; S' R] must be positive semidefinite ",
                 "(smallest eigenvalue of Q - S R^-1 S' is ",
                 signif(min_eig_sym(Sc), 3), "). Use ",
                 "`check_weights = FALSE` to solve an indefinite equation.")
    }
  } else if (rcond(R) < 1e-14) {
    ctrl_abort("`R` must be nonsingular.")
  }
  list(A = A, B = B, Q = Q, R = R, S = S, n = n, m = m)
}

riccati_solve <- function(A, B, Q, R, S, discrete, check_weights) {
  w <- riccati_inputs(A, B, Q, R, S, check_weights)
  n <- w$n
  if (n == 0L) return(riccati_core(w$A, w$B, w$Q, w$R, w$S, discrete))
  # Equilibrate before solving; the gain and the closed loop do not change.
  # 1. a common weight factor, which scales P only. It is applied first:
  #    dividing Q and R by rho multiplies G = B R^-1 B' by rho and divides Q
  #    by rho, which the state balancing below then evens out (applied
  #    last it would undo that balancing, for example after a change of
  #    input units)
  rho <- max(norm(w$Q, "1"), norm(w$R, "1"), norm(w$S, "1"))
  rho <- if (is.finite(rho) && rho > 0) 2^round(log2(rho)) else 1
  Q <- w$Q / rho
  R <- w$R / rho
  S <- w$S / rho
  # 2. symplectic diagonal state scaling x = T xb from balancing the
  #    Hamiltonian matrix, so that states in very different units carry
  #    comparable weight and G and Q have comparable size
  #    (Ab = T^-1 A T, Bb = T^-1 B, Qb = T Q T, Sb = T S)
  G <- w$B %*% solve(R, t(w$B))
  H <- rbind(cbind(w$A, -G), cbind(-Q, -t(w$A)))
  bs <- tryCatch(expm::balance(H, job = "S")$scale, error = function(e) rep(1, 2 * n))
  e <- round(log2(sqrt(bs[seq_len(n)] / bs[n + seq_len(n)])))
  e[!is.finite(e)] <- 0
  t <- 2^pmin(pmax(e, -100), 100)
  A <- w$A * outer(1 / t, t)
  B <- w$B / t
  Q <- Q * outer(t, t)
  S <- S * t
  # 3. continuous time: a common time scale (A, B and the weights divided
  #    by alpha leave P and K unchanged)
  alpha <- 1
  if (!discrete) {
    na <- norm(A, "1")
    if (is.finite(na) && na > 0) alpha <- 2^round(log2(na))
    A <- A / alpha
    B <- B / alpha
    Q <- Q / alpha
    R <- R / alpha
    S <- S / alpha
  }
  # 4. input units u = Du ub (B Du, Du R Du, S Du; P unchanged, K = Du Kb):
  #    R grows quadratically and B, S linearly, so du_j = |[B; S]_j| / R_jj
  #    balances each input column [B; -S; R] of the extended pencil. It
  #    comes last because it leaves G = B R^-1 B' and P unchanged, while
  #    the steps above rescale R and B differently.
  bsn <- sqrt(colSums(B^2) + colSums(S^2))
  rd <- diag(R)
  du <- ifelse(bsn > 0, bsn / rd, 1 / sqrt(rd))
  du <- 2^pmin(pmax(round(log2(du)), -200), 200)
  du[!is.finite(du)] <- 1
  B <- sweep(B, 2L, du, `*`)
  S <- sweep(S, 2L, du, `*`)
  R <- R * outer(du, du)
  out <- riccati_core(A, B, Q, R, S, discrete)
  out$P <- rho * out$P / outer(t, t)
  out$K <- du * out$K * rep(1 / t, each = nrow(out$K))
  out$poles <- out$poles * alpha
  out
}

riccati_core <- function(A, B, Q, R, S, discrete) {
  n <- nrow(A)
  m <- ncol(B)
  if (n == 0L) {
    return(structure(list(P = matrix(0, 0, 0), K = matrix(0, m, 0),
                          poles = complex(), residual = 0, symmetry_error = 0,
                          rcond_u11 = 1, discrete = discrete),
                     class = "ctrl_riccati"))
  }
  In <- diag(n)
  Z0 <- matrix(0, n, n)
  if (!discrete) {
    M <- rbind(cbind(A, Z0, B), cbind(-Q, -t(A), -S), cbind(t(S), t(B), R))
    N <- blkdiag(In, In, matrix(0, m, m))
  } else {
    M <- rbind(cbind(A, Z0, B), cbind(-Q, In, -S),
               cbind(t(S), matrix(0, m, n), R))
    N <- rbind(cbind(In, Z0, matrix(0, n, m)), cbind(Z0, t(A), matrix(0, n, m)),
               cbind(matrix(0, m, n), -t(B), matrix(0, m, m)))
  }
  # Left null space compression of the input block column [B; -S; R].
  X <- M[, 2L * n + seq_len(m), drop = FALSE]
  qx <- qr(X)
  if (qx$rank < m) {
    ctrl_abort("The input block of the Riccati pencil is rank deficient.")
  }
  W <- qr.Q(qx, complete = TRUE)[, m + seq_len(2L * n), drop = FALSE]
  Mr <- t(W) %*% M[, seq_len(2L * n), drop = FALSE]
  Nr <- t(W) %*% N[, seq_len(2L * n), drop = FALSE]
  # balance the pencil scale for the ordering test
  sc <- max(norm(Mr, "1"), norm(Nr, "1"), 1)
  tol_sel <- 100 * .Machine$double.eps * (2 * n)
  qz <- qz_real(Mr / sc, Nr / sc,
                order = if (discrete) "discrete" else "continuous",
                tol = tol_sel)
  if (qz$sdim != n) {
    ctrl_abort("No stabilizing solution: the ", if (discrete) "symplectic"
               else "Hamiltonian", " pencil has ", qz$sdim, " stable ",
               "eigenvalues instead of ", n, ". Eigenvalues on (or numerically ",
               "near) the ", if (discrete) "unit circle" else "imaginary axis",
               " usually mean that (A, B) is not stabilizable or that ",
               "(Q, A) has unobservable modes on the stability boundary.")
  }
  U11 <- qz$Z[seq_len(n), seq_len(n), drop = FALSE]
  U21 <- qz$Z[n + seq_len(n), seq_len(n), drop = FALSE]
  rc <- rcond(U11)
  if (rc < 1e-14) {
    ctrl_abort("No stabilizing solution: the stable deflating subspace is ",
               "not a graph subspace (rcond(U11) = ", signif(rc, 3), ").")
  }
  P <- t(solve(t(U11), t(U21)))
  sym_err <- asymmetry(P)
  P <- symmetrize(P)
  if (!discrete) {
    K <- solve(R, t(B) %*% P + t(S))
    Res <- t(A) %*% P + P %*% A - (P %*% B + S) %*% K + Q
    res <- rel_residual(Res, t(A) %*% P, Q, (P %*% B + S) %*% K)
  } else {
    G <- R + t(B) %*% P %*% B
    K <- solve(G, t(B) %*% P %*% A + t(S))
    Res <- t(A) %*% P %*% A - P - (t(A) %*% P %*% B + S) %*% K + Q
    res <- rel_residual(Res, t(A) %*% P %*% A, P, Q)
  }
  if (!is.finite(res) || res > 1e-6) {
    ctrl_abort("The computed Riccati solution has relative residual ",
               signif(res, 3), " (above 1e-6); the problem is too ill ",
               "conditioned for an accurate solution.")
  }
  poles <- as.complex(eigen(A - B %*% K, only.values = TRUE)$values)
  stable <- if (discrete) all(Mod(poles) < 1) else all(Re(poles) < 0)
  if (!stable) {
    ctrl_abort("The computed Riccati solution is not stabilizing (closed-loop ",
               "spectral ", if (discrete) "radius " else "abscissa ",
               signif(if (discrete) max(Mod(poles)) else max(Re(poles)), 4),
               "); the problem is numerically too close to the boundary.")
  }
  out <- list(P = P, K = K, poles = sort_poles(poles, as.numeric(discrete)),
              residual = res, symmetry_error = sym_err, rcond_u11 = rc,
              discrete = discrete)
  class(out) <- "ctrl_riccati"
  out
}

#' @export
print.ctrl_riccati <- function(x, digits = 4, ...) {
  cat(if (x$discrete) "Discrete" else "Continuous",
      " algebraic Riccati equation: stabilizing solution\n", sep = "")
  cat("P =\n")
  print(signif(x$P, digits))
  cat("relative residual: ", signif(x$residual, 3),
      ", symmetry error before symmetrization: ", signif(x$symmetry_error, 3),
      "\n", sep = "")
  cat("closed-loop poles: ", paste(format(signif(x$poles, digits)),
                                   collapse = ", "), "\n", sep = "")
  invisible(x)
}
