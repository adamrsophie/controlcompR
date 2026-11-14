#' Pole placement
#'
#' Compute a state-feedback gain `K` such that the eigenvalues of
#' \eqn{A - B K} are the requested `poles`.
#'
#' The method is the robust eigenstructure assignment of Kautsky, Nichols
#' and Van Dooren (method 0). With \eqn{B = U_0 Z} (QR), each closed-loop
#' eigenvector must lie in \eqn{S_i = \mathrm{null}(U_1^\top (A - \lambda_i
#' I))}; the eigenvectors are chosen iteratively to make the eigenvector
#' matrix `X` as well conditioned as possible, and
#' \eqn{K = Z^{-1} U_0^\top (A - X \Lambda X^{-1})}. For single-input
#' systems the eigenvectors are unique, so the result equals Ackermann's
#' formula, but it is computed without the controllability matrix.
#'
#' Requirements: `(A, B)` controllable, `B` of full column rank, complex
#' poles in conjugate pairs, and no pole repeated more than `rank(B)`
#' times. The achieved closed-loop poles are verified; `cond_X` (condition
#' number of the eigenvector matrix) measures how sensitive the placement is
#' to perturbations; large values mean the placed poles are fragile.
#'
#' @param sys A `ctrl_ss` model or the matrix `A`.
#' @param B Input matrix when `sys` is a matrix.
#' @param poles Desired closed-loop poles (length `n`).
#' @param maxiter Maximum number of KNV sweeps.
#' @return The `m x n` gain matrix `K` with attributes `poles` (achieved),
#'   `cond_X` and `max_error` (largest distance between requested and
#'   achieved poles).
#' @references Kautsky, J., Nichols, N. K. and Van Dooren, P. (1985). Robust
#'   pole assignment in linear state feedback. *International Journal of
#'   Control*, 41(5), 1129-1155.
#' @examples
#' plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
#' K <- ctrl_place(plant, poles = c(-4, -5))
#' eigen(plant$A - plant$B %*% K)$values
#' @export
ctrl_place <- function(sys, B = NULL, poles, maxiter = 50L) {
  # ctrl_place(model, poles): the second positional argument is the poles
  if (missing(poles) && is_ctrl_system(sys) && !is.null(B)) {
    poles <- B
    B <- NULL
  }
  ab <- get_ab(sys, B, "B")
  A <- ab$A
  B <- ab$B
  n <- nrow(A)
  m <- ncol(B)
  poles <- as.complex(poles)
  if (length(poles) != n) ctrl_abort("`poles` must have length ", n, ".")
  check_conjugate(poles, "poles")
  qb <- qr(B)
  if (qb$rank < m) ctrl_abort("`B` must have full column rank.")
  ctr <- ctrl_controllability(A, B)
  if (!ctr$full) {
    ctrl_abort("(A, B) is not controllable (rank ", ctr$rank, " of ", n,
               "); the uncontrollable modes cannot be moved.")
  }
  mult <- vapply(poles, function(p) sum(Mod(poles - p) < 1e-10 * max(1, Mod(p))),
                 0)
  if (any(mult > m)) {
    ctrl_abort("A pole is repeated more than rank(B) = ", m, " times.")
  }
  Qf <- qr.Q(qb, complete = TRUE)
  U0 <- Qf[, seq_len(m), drop = FALSE]
  U1 <- Qf[, seq_len(n - m) + m, drop = FALSE]
  Z <- qr.R(qb)[seq_len(m), seq_len(m), drop = FALSE]
  # order poles so that conjugate pairs are adjacent (positive imag first)
  ord <- order(Re(poles), -abs(Im(poles)), -Im(poles))
  poles <- poles[ord]
  real_pole <- abs(Im(poles)) <= 1e-12 * pmax(1, Mod(poles))
  poles[real_pole] <- Re(poles[real_pole])
  bases <- lapply(seq_len(n), function(i) {
    if (n == m) return(diag(n) + 0i)
    M <- t(U1) %*% (A - poles[i] * diag(n))
    if (real_pole[i]) M <- Re(M)
    s <- svd(M, nu = 0, nv = n)
    r <- sum(s$d > max(dim(M)) * .Machine$double.eps * max(s$d, 1))
    s$v[, seq.int(r + 1L, n), drop = FALSE]
  })
  partner <- rep(NA_integer_, n)
  i <- 1L
  while (i <= n) {
    if (!real_pole[i]) {
      partner[i] <- i + 1L
      partner[i + 1L] <- i
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  X <- matrix(0 + 0i, n, n)
  seen <- list()
  for (i in seq_len(n)) {
    key <- format(poles[i])
    k <- (seen[[key]] %||% 0L) + 1L
    seen[[key]] <- k
    S <- bases[[i]]
    X[, i] <- S[, min(k, ncol(S))]
  }
  for (i in seq_len(n)) if (!is.na(partner[i]) && partner[i] < i) X[, i] <- Conj(X[, partner[i]])
  cond_of <- function(X) {
    d <- svd(X, 0, 0)$d
    if (min(d) == 0) Inf else max(d) / min(d)
  }
  best <- cond_of(X)
  if (m > 1L) {
    for (iter in seq_len(maxiter)) {
      for (i in seq_len(n)) {
        if (!is.na(partner[i]) && partner[i] < i) next
        others <- X[, -i, drop = FALSE]
        qo <- qr(others)
        y <- qr.Q(qo, complete = TRUE)[, n]
        S <- bases[[i]]
        if (real_pole[i]) {
          S <- Re(S)
          v <- S %*% Re(t(S) %*% y)
          if (sqrt(sum(Mod(v)^2)) < 1e-14) v <- S %*% Im(t(S) %*% y)
        } else {
          v <- S %*% (Conj(t(S)) %*% y)
        }
        nv <- sqrt(sum(Mod(v)^2))
        if (nv > 1e-14) {
          X[, i] <- v / nv
          if (!is.na(partner[i])) X[, partner[i]] <- Conj(X[, i])
        }
      }
      cnew <- cond_of(X)
      if (!is.finite(cnew) || cnew >= best * (1 - 1e-6)) {
        best <- min(best, cnew)
        break
      }
      best <- cnew
    }
  }
  if (!is.finite(cond_of(X))) {
    ctrl_abort("Pole placement failed: the eigenvector matrix is singular.")
  }
  Lam <- diag(poles, n)
  M <- A - X %*% Lam %*% solve(X)
  K <- solve(Z, t(U0) %*% M)
  if (max(Mod(Im(K))) > 1e-6 * max(1, Mod(K))) {
    ctrl_abort("Pole placement produced a complex gain; check that complex ",
               "poles come in conjugate pairs.")
  }
  K <- Re(K)
  achieved <- as.complex(eigen(A - B %*% K, only.values = TRUE)$values)
  err <- max(vapply(poles, function(p) min(Mod(achieved - p)), 0))
  tol_err <- 1e-6 * max(1, Mod(poles), norm(A, "1"))
  if (err > tol_err) {
    ctrl_warn("Placed poles deviate from the requested poles by up to ",
              signif(err, 3), "; the placement is numerically sensitive ",
              "(cond(X) = ", signif(cond_of(X), 3), ").")
  }
  attr(K, "poles") <- achieved
  attr(K, "cond_X") <- cond_of(X)
  attr(K, "max_error") <- err
  K
}
