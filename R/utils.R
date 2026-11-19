#' @useDynLib controlcompR, .registration = TRUE
#' @importFrom stats rnorm setNames
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"

`%||%` <- function(x, y) if (is.null(x)) y else x

ctrl_abort <- function(...) stop(paste0(...), call. = FALSE)

ctrl_warn <- function(...) warning(paste0(...), call. = FALSE)

# Coerce a numeric input to a real matrix with explicit dimensions.
# `nrow`/`ncol` are only used to shape scalars and vectors, never to recycle.
as_real_matrix <- function(x, name, nrow = NULL, ncol = NULL) {
  if (is.null(x)) {
    if (is.null(nrow) || is.null(ncol)) ctrl_abort("`", name, "` is required.")
    return(matrix(0, nrow, ncol))
  }
  if (inherits(x, "Matrix")) x <- as.matrix(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.numeric(x) && !is.logical(x)) {
    ctrl_abort("`", name, "` must be numeric.")
  }
  if (is.complex(x)) ctrl_abort("`", name, "` must be real.")
  if (is.null(dim(x))) {
    if (length(x) == 1L && !is.null(nrow) && !is.null(ncol) &&
        nrow == 1L && ncol == 1L) {
      x <- matrix(x, 1, 1)
    } else if (length(x) == 1L && is.null(nrow) && is.null(ncol)) {
      x <- matrix(x, 1, 1)
    } else if (!is.null(nrow) && !is.null(ncol) && length(x) == nrow * ncol &&
               (nrow == 1L || ncol == 1L)) {
      x <- matrix(x, nrow, ncol)
    } else if (!is.null(ncol) && is.null(nrow) && length(x) == ncol) {
      x <- matrix(x, 1, ncol)
    } else if (!is.null(nrow) && is.null(ncol) && length(x) == nrow) {
      x <- matrix(x, nrow, 1)
    } else {
      ctrl_abort("`", name, "` must be a matrix; got a vector of length ",
                 length(x), ".")
    }
  }
  if (length(dim(x)) != 2L) ctrl_abort("`", name, "` must be a 2-d matrix.")
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    ctrl_abort("`", name, "` contains non-finite values (NA, NaN or Inf).")
  }
  dimnames(x) <- NULL
  x
}

check_dims <- function(x, name, nrow, ncol) {
  if (nrow(x) != nrow || ncol(x) != ncol) {
    ctrl_abort("`", name, "` must be ", nrow, " x ", ncol, "; got ",
               nrow(x), " x ", ncol(x), ".")
  }
  invisible(x)
}

check_square <- function(x, name) {
  if (nrow(x) != ncol(x)) {
    ctrl_abort("`", name, "` must be square; got ", nrow(x), " x ", ncol(x), ".")
  }
  invisible(x)
}

check_scalar_number <- function(x, name, positive = FALSE, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    ctrl_abort("`", name, "` must be a single finite number.")
  }
  if (positive && x <= 0) ctrl_abort("`", name, "` must be positive.")
  invisible(x)
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    ctrl_abort("`", name, "` must be TRUE or FALSE.")
  }
  invisible(x)
}

symmetrize <- function(X) (X + t(X)) / 2

# Relative asymmetry of a square matrix.
asymmetry <- function(X) {
  nx <- max(norm(X, "F"), .Machine$double.xmin)
  norm(X - t(X), "F") / nx
}

# Check that a matrix is symmetric (to a relative tolerance) and return the
# symmetrized matrix.
check_symmetric <- function(X, name, tol = 1e-10) {
  check_square(X, name)
  if (length(X) && asymmetry(X) > tol) {
    ctrl_abort("`", name, "` must be symmetric (relative asymmetry ",
               signif(asymmetry(X), 3), ").")
  }
  symmetrize(X)
}

# Smallest eigenvalue of a symmetric matrix, scaled by its norm.
min_eig_sym <- function(X) {
  if (!length(X)) return(Inf)
  min(eigen(symmetrize(X), symmetric = TRUE, only.values = TRUE)$values)
}

check_psd <- function(X, name, tol = NULL) {
  X <- check_symmetric(X, name)
  if (!length(X)) return(X)
  tol <- tol %||% (100 * nrow(X) * .Machine$double.eps * norm(X, "2"))
  if (min_eig_sym(X) < -tol) {
    ctrl_abort("`", name, "` must be positive semidefinite (smallest eigenvalue ",
               signif(min_eig_sym(X), 3), ").")
  }
  X
}

check_pd <- function(X, name, tol = NULL) {
  X <- check_symmetric(X, name)
  if (!length(X)) return(X)
  tol <- tol %||% (100 * nrow(X) * .Machine$double.eps * norm(X, "2"))
  if (min_eig_sym(X) <= tol) {
    ctrl_abort("`", name, "` must be positive definite (smallest eigenvalue ",
               signif(min_eig_sym(X), 3), ").")
  }
  X
}

# Scale-aware numerical rank via the SVD.
numerical_rank <- function(M, tol = NULL) {
  if (!length(M)) return(list(rank = 0L, sv = numeric(), tol = 0))
  sv <- svd(M, nu = 0, nv = 0)$d
  tol <- tol %||% (max(dim(M)) * .Machine$double.eps * max(sv, 0))
  list(rank = sum(sv > tol), sv = sv, tol = tol)
}

# Orthonormal basis of the null space of M (columns).
null_space <- function(M, tol = NULL) {
  n <- ncol(M)
  if (nrow(M) == 0L) return(diag(n))
  s <- svd(M, nu = 0, nv = n)
  tol <- tol %||% (max(dim(M)) * .Machine$double.eps * max(s$d, 0))
  r <- sum(s$d > tol)
  s$v[, seq.int(r + 1L, length.out = n - r), drop = FALSE]
}

# Characteristic polynomial coefficients (descending powers) from roots.
poly_from_roots <- function(r) {
  p <- 1 + 0i
  for (z in r) p <- c(p, 0) - z * c(0, p)
  if (all(abs(Im(p)) <= 1e-9 * pmax(1, abs(Re(p))))) Re(p) else p
}

# Roots of a polynomial given by descending coefficients.
poly_roots <- function(p) {
  p <- strip_leading_zeros(p)
  n <- length(p) - 1L
  if (n < 1L) return(complex())
  # trailing zeros are roots at the origin
  nz <- 0L
  while (n - nz >= 1L && p[length(p) - nz] == 0) nz <- nz + 1L
  q <- p[seq_len(length(p) - nz)]
  m <- length(q) - 1L
  r <- if (m >= 1L) {
    if (m == 1L) {
      -q[2] / q[1] + 0i
    } else {
      comp <- matrix(0, m, m)
      comp[1, ] <- -q[-1] / q[1]
      comp[cbind(2:m, 1:(m - 1))] <- 1
      as.complex(eigen(comp, only.values = TRUE)$values)
    }
  } else complex()
  c(r, rep(0 + 0i, nz))
}

strip_leading_zeros <- function(p, tol = 0) {
  if (!length(p)) return(0)
  scale <- max(abs(p), 0)
  keep <- which(abs(p) > tol * scale)
  if (!length(keep)) return(0)
  p[seq.int(keep[1], length(p))]
}

poly_eval <- function(p, z) {
  out <- rep(0 + 0i, length(z))
  for (a in p) out <- out * z + a
  out
}

poly_mul <- function(a, b) {
  if (!length(a) || !length(b)) return(0)
  out <- rep(0, length(a) + length(b) - 1L)
  for (i in seq_along(a)) {
    idx <- i:(i + length(b) - 1L)
    out[idx] <- out[idx] + a[i] * b
  }
  out
}

poly_add <- function(a, b) {
  n <- max(length(a), length(b))
  c(rep(0, n - length(a)), a) + c(rep(0, n - length(b)), b)
}

# Format a polynomial in descending powers of `var`.
poly_to_string <- function(p, var = "s", digits = 4) {
  p <- strip_leading_zeros(p)
  n <- length(p) - 1L
  terms <- character()
  for (i in seq_along(p)) {
    a <- p[i]
    if (a == 0) next
    k <- n - i + 1L
    mag <- abs(a)
    coef <- formatC(mag, digits = digits, format = "g")
    pw <- if (k == 0L) "" else if (k == 1L) var else paste0(var, "^", k)
    body <- if (k == 0L) coef else if (mag == 1) pw else paste0(coef, " ", pw)
    sgn <- if (a < 0) "-" else "+"
    terms <- c(terms, if (!length(terms)) paste0(if (a < 0) "-" else "", body)
                      else paste(sgn, body))
  }
  if (!length(terms)) "0" else paste(terms, collapse = " ")
}

# Default tolerance for eigenvalue classification relative to the matrix scale.
default_eig_tol <- function(A) {
  if (!length(A)) return(0)
  1e3 * nrow(A) * .Machine$double.eps * max(1, norm(A, "1"))
}

# Frobenius-norm relative residual.
rel_residual <- function(R, ...) {
  scale <- max(vapply(list(...), function(M) norm(as.matrix(M), "F"), 0), 1e-300)
  norm(R, "F") / scale
}

match_time_domain <- function(dt) {
  if (is.null(dt)) return(0)
  if (!is.numeric(dt) || length(dt) != 1L || is.na(dt) || !is.finite(dt) ||
      dt < 0) {
    ctrl_abort("`dt` must be NULL or 0 (continuous time) or a positive ",
               "finite sample interval (discrete time).")
  }
  as.numeric(dt)
}

default_names <- function(prefix, k, given = NULL, what = prefix) {
  if (is.null(given)) {
    if (k == 0L) return(character())
    return(paste0(prefix, seq_len(k)))
  }
  given <- as.character(given)
  if (length(given) != k) {
    ctrl_abort("`", what, "_names` must have length ", k, "; got ",
               length(given), ".")
  }
  if (anyNA(given) || any(!nzchar(given))) {
    ctrl_abort("`", what, "_names` must be non-empty strings.")
  }
  if (anyDuplicated(given)) ctrl_abort("`", what, "_names` must be unique.")
  given
}

blkdiag <- function(...) {
  mats <- list(...)
  if (length(mats) == 1L && is.list(mats[[1]]) && !is.matrix(mats[[1]])) {
    mats <- mats[[1]]
  }
  mats <- lapply(mats, as.matrix)
  nr <- vapply(mats, nrow, 0L)
  nc <- vapply(mats, ncol, 0L)
  out <- matrix(0, sum(nr), sum(nc))
  r0 <- 0L
  c0 <- 0L
  for (i in seq_along(mats)) {
    if (nr[i] && nc[i]) {
      out[r0 + seq_len(nr[i]), c0 + seq_len(nc[i])] <- mats[[i]]
    }
    r0 <- r0 + nr[i]
    c0 <- c0 + nc[i]
  }
  out
}

# Solve with an explicit conditioning check instead of a silent inverse.
safe_solve <- function(A, B, what = "matrix", rcond_tol = 1e-13) {
  if (!length(A)) return(B)
  rc <- rcond(A)
  if (!is.finite(rc) || rc < rcond_tol) {
    ctrl_abort("The ", what, " is singular or ill-conditioned (reciprocal ",
               "condition number ", signif(rc, 3), ").")
  }
  solve(A, B)
}

# Cholesky factor of a symmetric positive semidefinite matrix via the
# eigendecomposition (returns F with F %*% t(F) = X).
psd_sqrt_factor <- function(X) {
  if (!length(X)) return(matrix(0, 0, 0))
  e <- eigen(symmetrize(X), symmetric = TRUE)
  vals <- pmax(e$values, 0)
  keep <- vals > max(vals, 0) * nrow(X) * .Machine$double.eps
  e$vectors[, keep, drop = FALSE] %*% diag(sqrt(vals[keep]), sum(keep))
}

# OSQP prints a polishing note even with verbose = FALSE; keep it quiet.
osqp_solve <- function(model) {
  res <- NULL
  utils::capture.output(res <- model@Solve())
  res
}

# Diagonal state-coordinate scaling that balances [A B; C 0] (LAPACK dgebal,
# scaling only). The transfer function is unchanged; results that depend on
# conditioning (frequency response, Hamiltonian tests) no longer depend on
# the units of the states.
prescale_ss <- function(s) {
  n <- nrow(s$A)
  if (!n) return(s)
  m <- ncol(s$B)
  p <- nrow(s$C)
  k <- max(m, p)
  M <- matrix(0, n + k, n + k)
  M[seq_len(n), seq_len(n)] <- s$A
  M[seq_len(n), n + seq_len(m)] <- s$B
  M[n + seq_len(p), seq_len(n)] <- s$C
  if (!all(is.finite(M))) return(s)
  d <- expm::balance(M, job = "S")$scale[seq_len(n)]
  if (!all(is.finite(d)) || any(d <= 0)) return(s)
  s$A <- s$A * outer(1 / d, d)
  s$B <- s$B / d
  s$C <- sweep(s$C, 2L, d, `*`)
  s
}

# Measurement matrix for the estimators: numeric, one column per output,
# rows complete or entirely NA (missing); infinite or NaN values are errors.
measurement_matrix <- function(y, p) {
  Y <- if (is.null(dim(y))) matrix(y, ncol = 1) else as.matrix(y)
  if (!(is.numeric(Y) || (is.logical(Y) && all(is.na(Y))))) {
    ctrl_abort("`y` must be numeric (NA marks missing samples).")
  }
  storage.mode(Y) <- "double"
  if (ncol(Y) != p) ctrl_abort("`y` must have ", p, " columns.")
  if (any(is.infinite(Y)) || any(is.nan(Y))) {
    ctrl_abort("`y` contains infinite or NaN values; only NA marks missing data.")
  }
  part <- apply(Y, 1, function(r) any(is.na(r)) && !all(is.na(r)))
  if (any(part)) {
    ctrl_abort("Rows of `y` must be complete or entirely missing (NA); ",
               "partially missing rows are not supported.")
  }
  Y
}
