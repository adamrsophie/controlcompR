#' Transfer-function models
#'
#' Construct a rational transfer function
#' \deqn{G(s) = \frac{b_m s^m + \cdots + b_1 s + b_0}{a_n s^n + \cdots +
#' a_1 s + a_0}}
#' in continuous time (`dt = NULL` or `0`) or in the variable `z` for a
#' discrete time model with sample interval `dt`.
#'
#' **Coefficient ordering.** Polynomials are given in *descending* powers,
#' as in MATLAB and python-control: `c(1, 3, 2)` is \eqn{s^2 + 3s + 2}.
#' Leading zeros are removed and the denominator is normalized to be monic;
#' the numerator is scaled by the same factor.
#'
#' **MIMO models.** For `p` outputs and `m` inputs, pass `num` as a `p x m`
#' list matrix of numerator vectors and `den` either as a single vector
#' (common denominator) or as a `p x m` list matrix.
#'
#' **Improper models.** A transfer function with a numerator degree above
#' the denominator degree can be created (for example `ctrl_tf("s")`, which
#' is handy for building models with arithmetic) and evaluated in the
#' frequency domain, but it has no finite-dimensional state-space
#' realization and [ctrl_tf2ss()] refuses to convert it.
#'
#' @param num Numerator coefficients (descending powers), a list matrix for
#'   MIMO models, or the string `"s"` or `"z"` for the Laplace or shift
#'   variable.
#' @param den Denominator coefficients (descending powers) or a list matrix.
#' @inheritParams ctrl_ss
#' @return An object of class `c("ctrl_tf", "ctrl_system")` with fields
#'   `num` and `den` (list matrices), `dt` and signal names.
#' @seealso [ctrl_zpk()], [ctrl_tf2ss()], [ctrl_ss2tf()]
#' @examples
#' G <- ctrl_tf(1, c(1, 3, 2))
#' G
#' s <- ctrl_tf("s")
#' H <- (s + 1) / (s^2 + 2 * s + 5)
#' H
#' # A discrete time integrator with sample interval 0.1
#' ctrl_tf(c(0.1, 0), c(1, -1), dt = 0.1)
#' @export
ctrl_tf <- function(num, den = 1, dt = NULL, input_names = NULL,
                    output_names = NULL) {
  dt <- match_time_domain(dt)
  if (is.character(num)) {
    if (length(num) != 1L || !num %in% c("s", "z")) {
      ctrl_abort("`num` must be numeric, a list matrix, or \"s\" / \"z\".")
    }
    if (num == "z" && dt == 0) {
      ctrl_abort("The shift variable \"z\" needs a positive `dt`.")
    }
    if (num == "s" && dt > 0) {
      ctrl_abort("The Laplace variable \"s\" is continuous time; use \"z\".")
    }
    num <- c(1, 0)
    den <- 1
  }
  if (is.list(num)) {
    if (is.null(dim(num))) num <- matrix(num, 1, length(num))
    p <- nrow(num)
    m <- ncol(num)
    if (!is.list(den)) {
      den <- matrix(rep(list(den), p * m), p, m)
    } else if (is.null(dim(den))) {
      den <- matrix(den, 1, length(den))
    }
    if (!identical(dim(den), dim(num))) {
      ctrl_abort("`den` must have the same dimensions as `num`.")
    }
  } else {
    num <- matrix(list(num), 1, 1)
    den <- matrix(list(den), 1, 1)
    p <- 1L
    m <- 1L
  }
  for (k in seq_along(num)) {
    nd <- normalize_tf(num[[k]], den[[k]])
    num[[k]] <- nd$num
    den[[k]] <- nd$den
  }
  structure(
    list(num = num, den = den, dt = dt,
         input_names = default_names("u", m, input_names, "input"),
         output_names = default_names("y", p, output_names, "output")),
    class = c("ctrl_tf", "ctrl_system")
  )
}

normalize_tf <- function(num, den) {
  check_poly(num, "num")
  check_poly(den, "den")
  den <- strip_leading_zeros(as.numeric(den))
  if (all(den == 0)) ctrl_abort("A denominator polynomial is identically zero.")
  num <- strip_leading_zeros(as.numeric(num))
  lead <- den[1]
  num <- num / lead
  den <- den / lead
  if (any(!is.finite(num)) || any(!is.finite(den))) {
    ctrl_abort("Making the denominator monic overflows double precision ",
               "(leading denominator coefficient ", signif(lead, 3), "). ",
               "Rescale the numerator and denominator.")
  }
  list(num = num, den = den)
}

check_poly <- function(p, name) {
  if (!is.numeric(p) || !length(p) || any(!is.finite(p))) {
    ctrl_abort("`", name, "` must be a non-empty finite numeric vector.")
  }
}

unique_denominators <- function(sys) {
  dens <- sys$den
  keys <- vapply(dens, function(d) paste(format(d, digits = 15), collapse = ","),
                 "")
  dens[!duplicated(keys)]
}

tf_degree <- function(p) length(strip_leading_zeros(p)) - 1L

#' Test whether a transfer function is proper
#'
#' @param sys A `ctrl_tf` model.
#' @param strict If `TRUE`, test strict properness.
#' @return `TRUE` or `FALSE`.
#' @examples
#' ctrl_is_proper(ctrl_tf(c(1, 0), c(1, 1)))
#' ctrl_is_proper(ctrl_tf(c(1, 0), c(1, 1)), strict = TRUE)
#' @export
ctrl_is_proper <- function(sys, strict = FALSE) {
  if (is_ctrl_ss(sys)) return(if (strict) all(sys$D == 0) else TRUE)
  if (!is_ctrl_tf(sys)) sys <- as_tf(sys)
  ok <- vapply(seq_along(sys$num), function(k) {
    dn <- tf_degree(sys$num[[k]])
    dd <- tf_degree(sys$den[[k]])
    if (all(sys$num[[k]] == 0)) TRUE else if (strict) dn < dd else dn <= dd
  }, TRUE)
  all(ok)
}

#' @export
print.ctrl_tf <- function(x, digits = 4, ...) {
  var <- if (x$dt > 0) "z" else "s"
  cat("Transfer function (", time_domain_label(x$dt), ")\n", sep = "")
  p <- nrow(x$num)
  m <- ncol(x$num)
  for (i in seq_len(p)) for (j in seq_len(m)) {
    if (p > 1L || m > 1L) {
      cat("\nFrom ", x$input_names[j], " to ", x$output_names[i], ":\n", sep = "")
    } else {
      cat("\n")
    }
    ns <- poly_to_string(x$num[[i, j]], var, digits)
    ds <- poly_to_string(x$den[[i, j]], var, digits)
    if (ds == "1") {
      cat("  ", ns, "\n", sep = "")
    } else {
      w <- max(nchar(ns), nchar(ds))
      pad <- function(s) paste0(strrep(" ", (w - nchar(s)) %/% 2), s)
      cat("  ", pad(ns), "\n  ", strrep("-", w), "\n  ", pad(ds), "\n", sep = "")
    }
  }
  invisible(x)
}

#' @export
summary.ctrl_tf <- function(object, ...) {
  out <- list(
    dims = ctrl_dims(object),
    dt = object$dt,
    poles = ctrl_poles(object),
    stability = ctrl_stability(object)$classification,
    controllability_rank = NA,
    observability_rank = NA,
    dc_gain = tryCatch(ctrl_dcgain(object), error = function(e) NULL)
  )
  class(out) <- "ctrl_summary"
  out
}

#' Zero-pole-gain models
#'
#' A SISO model \eqn{G(s) = k \prod (s - z_i) / \prod (s - p_i)}. Complex
#' zeros and poles must appear in conjugate pairs so that the model has
#' real coefficients.
#'
#' @param zeros Vector of zeros (possibly complex).
#' @param poles Vector of poles (possibly complex).
#' @param gain Real gain `k`.
#' @inheritParams ctrl_ss
#' @return An object of class `c("ctrl_zpk", "ctrl_system")`.
#' @examples
#' H <- ctrl_zpk(-1, c(-2, -3 + 1i, -3 - 1i), 4)
#' H
#' ctrl_zpk2tf(H)
#' @export
ctrl_zpk <- function(zeros = numeric(), poles = numeric(), gain = 1,
                     dt = NULL) {
  dt <- match_time_domain(dt)
  check_scalar_number(gain, "gain")
  zeros <- as.complex(zeros)
  poles <- as.complex(poles)
  check_conjugate(zeros, "zeros")
  check_conjugate(poles, "poles")
  structure(list(zeros = zeros, poles = poles, gain = as.numeric(gain),
                 dt = dt, input_names = "u1", output_names = "y1"),
            class = c("ctrl_zpk", "ctrl_system"))
}

check_conjugate <- function(z, name) {
  if (!length(z)) return(invisible())
  if (any(!is.finite(z))) ctrl_abort("`", name, "` must be finite.")
  tol <- 1e-8 * max(1, Mod(z))
  # greedy nearest matching of each value with the conjugate of another
  left <- seq_along(z)
  while (length(left)) {
    i <- left[1]
    if (abs(Im(z[i])) <= tol) {
      left <- left[-1]
      next
    }
    cand <- left[-1]
    d <- Mod(z[cand] - Conj(z[i]))
    if (!length(cand) || min(d) > tol) {
      ctrl_abort("Complex `", name, "` must appear in conjugate pairs.")
    }
    left <- setdiff(left, c(i, cand[which.min(d)]))
  }
}

sort_complex <- function(z) z[order(Re(z), Im(z))]

#' @export
print.ctrl_zpk <- function(x, digits = 4, ...) {
  cat("Zero-pole-gain model (", time_domain_label(x$dt), ")\n", sep = "")
  cat("  gain: ", signif(x$gain, digits), "\n", sep = "")
  cat("  zeros: ", if (length(x$zeros)) paste(format(signif(x$zeros, digits)),
                                              collapse = ", ") else "none",
      "\n", sep = "")
  cat("  poles: ", if (length(x$poles)) paste(format(signif(x$poles, digits)),
                                              collapse = ", ") else "none",
      "\n", sep = "")
  invisible(x)
}

#' Conversions between model representations
#'
#' * `ctrl_zpk2tf()` and `ctrl_tf2zpk()` convert between zero-pole-gain
#'   and polynomial forms (SISO).
#' * `ctrl_tf2ss()` realizes a proper transfer function in state space.
#'   SISO models use the controllable canonical form. MIMO models are
#'   realized element by element and then reduced with [ctrl_minreal()]
#'   when `minimal = TRUE`.
#' * `ctrl_ss2tf()` computes each SISO channel's transfer function. With
#'   `minimal = TRUE` (default) each channel is first reduced to a minimal
#'   realization, so uncontrollable or unobservable modes do not appear as
#'   canceling pole-zero pairs. The numerator of channel \eqn{(i,j)} is
#'   computed from the identity
#'   \eqn{c (sI - A)^{-1} b = [\det(sI - A + bc) - \det(sI - A)] /
#'   \det(sI - A)}.
#'
#' Root finding and characteristic polynomials are ill-conditioned for
#' high-order models; prefer state-space computations for large systems.
#'
#' @param sys A model.
#' @param minimal Reduce to a minimal realization (see details).
#' @param tol Tolerance passed to [ctrl_minreal()].
#' @return The converted model.
#' @examples
#' G <- ctrl_tf(c(1, 2), c(1, 3, 2))
#' ss <- ctrl_tf2ss(G)
#' ss
#' ctrl_ss2tf(ss)
#' ctrl_tf2zpk(G)
#' @export
ctrl_tf2ss <- function(sys, minimal = TRUE, tol = NULL) {
  if (is_ctrl_ss(sys)) return(sys)
  if (inherits(sys, "ctrl_zpk")) sys <- ctrl_zpk2tf(sys)
  if (!is_ctrl_tf(sys)) ctrl_abort("`sys` must be a transfer function.")
  if (!ctrl_is_proper(sys)) {
    ctrl_abort("The transfer function is improper (numerator degree above ",
               "denominator degree) and has no finite-dimensional state-space ",
               "realization. Add poles (for example a derivative filter) to ",
               "make it proper.")
  }
  p <- nrow(sys$num)
  m <- ncol(sys$num)
  if (p == 1L && m == 1L) {
    out <- siso_tf2ss(sys$num[[1]], sys$den[[1]], sys$dt)
    out <- ctrl_ss(out$A, out$B, out$C, out$D, dt = sys$dt,
                   input_names = sys$input_names,
                   output_names = sys$output_names)
    # the controllable canonical form is minimal unless numerator and
    # denominator share a factor (then it is unobservable)
    if (minimal && nrow(out$A) && !ctrl_is_observable(out)) {
      out <- ctrl_minreal(out, tol = tol)
    }
    return(out)
  }
  parts <- list()
  for (i in seq_len(p)) for (j in seq_len(m)) {
    r <- siso_tf2ss(sys$num[[i, j]], sys$den[[i, j]], sys$dt)
    Bj <- matrix(0, nrow(r$A), m)
    Bj[, j] <- r$B
    Ci <- matrix(0, p, nrow(r$A))
    Ci[i, ] <- r$C
    Di <- matrix(0, p, m)
    Di[i, j] <- r$D
    parts[[length(parts) + 1L]] <- list(A = r$A, B = Bj, C = Ci, D = Di)
  }
  A <- blkdiag(lapply(parts, `[[`, "A"))
  B <- do.call(rbind, lapply(parts, `[[`, "B"))
  C <- do.call(cbind, lapply(parts, `[[`, "C"))
  D <- Reduce(`+`, lapply(parts, `[[`, "D"))
  out <- if (nrow(A)) {
    ctrl_ss(A, B, C, D, dt = sys$dt, input_names = sys$input_names,
            output_names = sys$output_names)
  } else {
    ctrl_gain(D, dt = sys$dt, input_names = sys$input_names,
              output_names = sys$output_names)
  }
  if (minimal && nrow(out$A)) out <- ctrl_minreal(out, tol = tol)
  out
}

# Controllable canonical realization of a proper SISO transfer function.
siso_tf2ss <- function(num, den, dt) {
  den <- strip_leading_zeros(den)
  num <- strip_leading_zeros(num)
  n <- length(den) - 1L
  den <- den / den[1]
  if (all(num == 0)) {
    return(list(A = matrix(0, 0, 0), B = matrix(0, 0, 1), C = matrix(0, 1, 0),
                D = matrix(0, 1, 1)))
  }
  numf <- c(rep(0, n + 1L - length(num)), num)
  d <- numf[1]
  if (n == 0L) {
    return(list(A = matrix(0, 0, 0), B = matrix(0, 0, 1), C = matrix(0, 1, 0),
                D = matrix(d, 1, 1)))
  }
  # strictly proper remainder: num - d * den
  rem <- numf - d * den
  A <- matrix(0, n, n)
  A[1, ] <- -den[-1]
  if (n > 1L) A[cbind(2:n, 1:(n - 1L))] <- 1
  B <- matrix(c(1, rep(0, n - 1L)), n, 1)
  C <- matrix(rem[-1], 1, n)
  list(A = A, B = B, C = C, D = matrix(d, 1, 1))
}

#' @rdname ctrl_tf2ss
#' @export
ctrl_ss2tf <- function(sys, minimal = TRUE, tol = NULL) {
  if (is_ctrl_tf(sys)) return(sys)
  if (inherits(sys, "ctrl_zpk")) return(ctrl_zpk2tf(sys))
  sys <- as_ss(sys)
  p <- nrow(sys$C)
  m <- ncol(sys$B)
  num <- matrix(list(), p, m)
  den <- matrix(list(), p, m)
  for (i in seq_len(p)) for (j in seq_len(m)) {
    sub <- sys[i, j]
    if (minimal && nrow(sub$A)) sub <- ctrl_minreal(sub, tol = tol)
    nd <- siso_ss2tf(sub)
    num[[i, j]] <- nd$num
    den[[i, j]] <- nd$den
  }
  ctrl_tf(num, den, dt = sys$dt, input_names = sys$input_names,
          output_names = sys$output_names)
}

siso_ss2tf <- function(s) {
  n <- nrow(s$A)
  d <- s$D[1, 1]
  if (n == 0L) return(list(num = d, den = 1))
  charp <- function(M) Re(poly_from_roots(eigen(M, only.values = TRUE)$values))
  den <- charp(s$A)
  # C adj(sI - A) B = det(sI - A + alpha B C) / alpha - det(sI - A) / alpha
  # for any alpha; choosing alpha so that alpha B C is as large as A keeps
  # the difference well above rounding level for weak but genuine gains
  BC <- s$B %*% s$C
  g <- norm(BC, "1")
  # no input to output path through the states: only the feedthrough
  # remains (kept over the full denominator, as requested)
  if (g == 0) return(list(num = d * den, den = den))
  a <- norm(s$A, "1")
  alpha <- (if (a > 0) a else 1) / g
  num <- charp(s$A - alpha * BC) - den
  num[1] <- 0
  # remove leading coefficients that are pure rounding noise. Compare in
  # the frequency variable s / c, with c the size of A, so that a slow model
  # (small c) keeps its small low-order coefficients; logs avoid overflow
  cl <- log(if (a > 0) a else 1)
  pw <- (length(num) - 1):0
  lnum <- log(abs(num)) + (pw - n) * cl
  lscale <- max(lnum, log(abs(den)) + (pw - n) * cl)
  for (k in seq_along(num)) {
    if (lnum[k] <= log(1e3 * n * .Machine$double.eps) + lscale) num[k] <- 0 else break
  }
  num <- num / alpha + d * den
  list(num = strip_leading_zeros(num), den = den)
}

#' @rdname ctrl_tf2ss
#' @export
ctrl_zpk2tf <- function(sys) {
  if (!inherits(sys, "ctrl_zpk")) ctrl_abort("`sys` must be a `ctrl_zpk` model.")
  num <- Re(poly_from_roots(sys$zeros)) * sys$gain
  den <- Re(poly_from_roots(sys$poles))
  ctrl_tf(num, den, dt = sys$dt)
}

#' @rdname ctrl_tf2ss
#' @export
ctrl_tf2zpk <- function(sys) {
  if (!is_ctrl_tf(sys)) sys <- ctrl_ss2tf(sys)
  if (nrow(sys$num) != 1L || ncol(sys$num) != 1L) {
    ctrl_abort("Zero-pole-gain form is only available for SISO models.")
  }
  num <- sys$num[[1]]
  den <- sys$den[[1]]
  k <- num[1] / den[1]
  if (all(num == 0)) k <- 0
  ctrl_zpk(poly_roots(num), poly_roots(den), k, dt = sys$dt)
}

as_tf <- function(sys) {
  if (is_ctrl_tf(sys)) return(sys)
  if (inherits(sys, "ctrl_zpk")) return(ctrl_zpk2tf(sys))
  if (is_ctrl_ss(sys)) return(ctrl_ss2tf(sys))
  if (is.numeric(sys)) return(ctrl_tf(sys, 1))
  ctrl_abort("Cannot convert object to a transfer function.")
}

# Evaluate a transfer function at complex points; returns p x m x length(z).
tf_eval <- function(sys, z) {
  p <- nrow(sys$num)
  m <- ncol(sys$num)
  out <- array(0 + 0i, c(p, m, length(z)))
  for (i in seq_len(p)) for (j in seq_len(m)) {
    out[i, j, ] <- poly_eval(sys$num[[i, j]], z) / poly_eval(sys$den[[i, j]], z)
  }
  out
}

#' Pole-zero cancellation for transfer functions
#'
#' Remove numerator and denominator roots that coincide within `tol`
#' (relative to the root magnitude). Cancellations are a numerical
#' decision: a near cancellation of an unstable pole hides an internal
#' instability, so inspect the result before using it for design.
#'
#' @param sys A SISO `ctrl_tf` model.
#' @param tol Relative tolerance for matching roots.
#' @return A `ctrl_tf` model.
#' @examples
#' G <- ctrl_tf(c(1, 1), c(1, 3, 2))  # (s + 1) / ((s + 1)(s + 2))
#' ctrl_tf_minreal(G)
#' @export
ctrl_tf_minreal <- function(sys, tol = sqrt(.Machine$double.eps)) {
  if (!is_ctrl_tf(sys) || length(sys$num) != 1L) {
    ctrl_abort("`sys` must be a SISO `ctrl_tf` model.")
  }
  zpk <- ctrl_tf2zpk(sys)
  z <- zpk$zeros
  p <- zpk$poles
  keepz <- rep(TRUE, length(z))
  keepp <- rep(TRUE, length(p))
  for (i in seq_along(z)) {
    d <- Mod(p - z[i])
    d[!keepp] <- Inf
    j <- which.min(d)
    if (length(j) && d[j] <= tol * max(1, Mod(z[i]))) {
      keepz[i] <- FALSE
      keepp[j] <- FALSE
    }
  }
  ctrl_zpk2tf(ctrl_zpk(z[keepz], p[keepp], zpk$gain, dt = sys$dt))
}
