#' State-space models
#'
#' Construct a linear time-invariant state-space model
#' \deqn{\dot x(t) = A x(t) + B u(t), \quad y(t) = C x(t) + D u(t)}
#' in continuous time, or
#' \deqn{x_{k+1} = A x_k + B u_k, \quad y_k = C x_k + D u_k}
#' in discrete time with sample interval `dt`.
#'
#' The time domain is explicit: `dt = NULL` or `dt = 0` gives a continuous
#' time model, and a positive `dt` gives a discrete time model sampled every
#' `dt` time units. There is no "unknown sample time"; discrete models
#' always carry their sample interval.
#'
#' Matrices are validated for finiteness and for consistent dimensions:
#' `A` is `n x n`, `B` is `n x m`, `C` is `p x n` and `D` is `p x m`.
#' `D` defaults to a zero matrix. A model with `n = 0` is a static gain
#' `y = D u`; pass `A = NULL` together with `D`, or use [ctrl_gain()].
#' Scalars are accepted for `1 x 1` matrices and vectors for single-row or
#' single-column matrices only when the shape is unambiguous.
#'
#' @param A State matrix (`n x n`).
#' @param B Input matrix (`n x m`).
#' @param C Output matrix (`p x n`).
#' @param D Feedthrough matrix (`p x m`); defaults to zero.
#' @param dt `NULL` or `0` for continuous time, or a positive sample interval.
#' @param state_names Optional unique state names.
#' @param input_names Optional unique input names.
#' @param output_names Optional unique output names.
#'   Defaults are `x1, x2, ...`, `u1, ...` and `y1, ...`.
#' @return An object of class `c("ctrl_ss", "ctrl_system")`.
#' @seealso [ctrl_tf()], [ctrl_ssdata()], [ctrl_c2d()], [ctrl_poles()]
#' @examples
#' A <- matrix(c(0, 1, -2, -3), 2, byrow = TRUE)
#' plant <- ctrl_ss(A, B = c(0, 1), C = c(1, 0))
#' plant
#' ctrl_poles(plant)
#'
#' # A discrete time model sampled every 0.1 time units
#' ctrl_ss(0.9, 0.1, 1, dt = 0.1)
#' @export
ctrl_ss <- function(A, B = NULL, C = NULL, D = NULL, dt = NULL,
                    state_names = NULL, input_names = NULL,
                    output_names = NULL) {
  dt <- match_time_domain(dt)
  if (is.null(A) || length(A) == 0L) {
    if (is.null(D)) ctrl_abort("A static gain model (no states) needs `D`.")
    D <- as_real_matrix(D, "D")
    n <- 0L
    m <- ncol(D)
    p <- nrow(D)
    A <- matrix(0, 0, 0)
    B <- matrix(0, 0, m)
    C <- matrix(0, p, 0)
  } else {
    A <- as_real_matrix(A, "A")
    check_square(A, "A")
    n <- nrow(A)
    if (is.null(B)) ctrl_abort("`B` is required.")
    B <- as_real_matrix(B, "B", nrow = n)
    if (nrow(B) != n) {
      ctrl_abort("`B` must have ", n, " rows to match `A`; got ", nrow(B), ".")
    }
    m <- ncol(B)
    if (is.null(C)) ctrl_abort("`C` is required.")
    C <- as_real_matrix(C, "C", ncol = n)
    if (ncol(C) != n) {
      ctrl_abort("`C` must have ", n, " columns to match `A`; got ", ncol(C), ".")
    }
    p <- nrow(C)
    D <- if (is.null(D)) matrix(0, p, m) else as_real_matrix(D, "D", p, m)
    check_dims(D, "D", p, m)
  }
  sys <- structure(
    list(A = A, B = B, C = C, D = D, dt = dt,
         state_names = default_names("x", n, state_names, "state"),
         input_names = default_names("u", m, input_names, "input"),
         output_names = default_names("y", p, output_names, "output")),
    class = c("ctrl_ss", "ctrl_system")
  )
  sys
}

#' Static gain model
#'
#' A state-space model without states, `y = D u`.
#'
#' @param D Gain matrix (`p x m`) or scalar.
#' @inheritParams ctrl_ss
#' @return A `ctrl_ss` object with zero states.
#' @examples
#' ctrl_gain(diag(2))
#' @export
ctrl_gain <- function(D, dt = NULL, input_names = NULL, output_names = NULL) {
  ctrl_ss(NULL, D = D, dt = dt, input_names = input_names,
          output_names = output_names)
}

#' Test for control system objects
#'
#' @param x Any object.
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_ctrl_ss(ctrl_ss(-1, 1, 1))
#' is_ctrl_tf(ctrl_tf(1, c(1, 1)))
#' @export
is_ctrl_ss <- function(x) inherits(x, "ctrl_ss")

#' @rdname is_ctrl_ss
#' @export
is_ctrl_tf <- function(x) inherits(x, "ctrl_tf")

#' @rdname is_ctrl_ss
#' @export
is_ctrl_system <- function(x) inherits(x, "ctrl_system")

assert_system <- function(sys, name = "sys") {
  if (!is_ctrl_system(sys)) {
    ctrl_abort("`", name, "` must be a `ctrl_ss` or `ctrl_tf` model.")
  }
  invisible(sys)
}

# Coerce any supported model to state space.
as_ss <- function(sys, name = "sys") {
  if (is_ctrl_ss(sys)) return(sys)
  if (is_ctrl_tf(sys)) return(ctrl_tf2ss(sys))
  if (inherits(sys, "ctrl_zpk")) return(ctrl_tf2ss(ctrl_zpk2tf(sys)))
  ctrl_abort("`", name, "` must be a `ctrl_ss`, `ctrl_tf` or `ctrl_zpk` model.")
}

#' State-space data and dimensions
#'
#' Accessors that avoid depending on the internal fields of model objects.
#'
#' @param sys A `ctrl_ss`, `ctrl_tf` or `ctrl_zpk` model (transfer functions
#'   are converted with [ctrl_tf2ss()]).
#' @return `ctrl_ssdata()` returns a list with `A`, `B`, `C`, `D` and `dt`.
#'   `ctrl_dims()` returns a named integer vector with `states`, `inputs` and
#'   `outputs`. `ctrl_is_discrete()` and `ctrl_is_continuous()` return a
#'   logical, and `ctrl_sample_time()` returns `0` for continuous models.
#' @examples
#' sys <- ctrl_ss(-1, 1, 1)
#' ctrl_dims(sys)
#' ctrl_ssdata(sys)$A
#' ctrl_is_discrete(ctrl_c2d(sys, 0.1))
#' @export
ctrl_ssdata <- function(sys) {
  s <- as_ss(sys)
  list(A = s$A, B = s$B, C = s$C, D = s$D, dt = s$dt)
}

#' @rdname ctrl_ssdata
#' @export
ctrl_dims <- function(sys) {
  if (is_ctrl_tf(sys)) {
    return(c(states = NA_integer_, inputs = ncol(sys$num),
             outputs = nrow(sys$num)))
  }
  s <- as_ss(sys)
  c(states = nrow(s$A), inputs = ncol(s$B), outputs = nrow(s$C))
}

#' @rdname ctrl_ssdata
#' @export
ctrl_is_discrete <- function(sys) {
  assert_system(sys)
  sys$dt > 0
}

#' @rdname ctrl_ssdata
#' @export
ctrl_is_continuous <- function(sys) {
  assert_system(sys)
  sys$dt == 0
}

#' @rdname ctrl_ssdata
#' @export
ctrl_sample_time <- function(sys) {
  assert_system(sys)
  sys$dt
}

#' Signal names of a model
#'
#' @param sys A model.
#' @return A list with `states`, `inputs` and `outputs` character vectors.
#' @examples
#' ctrl_signal_names(ctrl_ss(-1, 1, 1, input_names = "force"))
#' @export
ctrl_signal_names <- function(sys) {
  assert_system(sys)
  list(states = sys$state_names %||% character(), inputs = sys$input_names,
       outputs = sys$output_names)
}

time_domain_label <- function(dt) {
  if (dt > 0) paste0("discrete time, dt = ", format(dt)) else "continuous time"
}

#' @export
print.ctrl_ss <- function(x, digits = 4, ...) {
  d <- ctrl_dims(x)
  cat("State-space model (", time_domain_label(x$dt), ")\n", sep = "")
  cat(d["states"], " states, ", d["inputs"], " inputs, ", d["outputs"],
      " outputs\n", sep = "")
  show_mat <- function(M, rn, cn, label) {
    cat("\n", label, " =\n", sep = "")
    if (!length(M)) {
      cat("  (empty ", nrow(M), " x ", ncol(M), ")\n", sep = "")
      return(invisible())
    }
    dimnames(M) <- list(rn, cn)
    print(signif(M, digits))
  }
  show_mat(x$A, x$state_names, x$state_names, "A")
  show_mat(x$B, x$state_names, x$input_names, "B")
  show_mat(x$C, x$output_names, x$state_names, "C")
  show_mat(x$D, x$output_names, x$input_names, "D")
  invisible(x)
}

#' Summaries of control models
#'
#' @param object A model.
#' @param ... Unused.
#' @return A `ctrl_summary` list with dimensions, time domain, poles,
#'   stability, controllability and observability ranks.
#' @examples
#' summary(ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#' @export
summary.ctrl_ss <- function(object, ...) {
  st <- ctrl_stability(object)
  n <- nrow(object$A)
  out <- list(
    dims = ctrl_dims(object),
    dt = object$dt,
    poles = st$poles,
    stability = st$classification,
    controllability_rank = if (n) ctrb_check(as_ss(object)$A, as_ss(object)$B)$rank else 0L,
    observability_rank = if (n) obsv_check(as_ss(object)$A, as_ss(object)$C)$rank else 0L,
    dc_gain = tryCatch(ctrl_dcgain(object), error = function(e) NULL)
  )
  class(out) <- "ctrl_summary"
  out
}

#' @export
print.ctrl_summary <- function(x, digits = 4, ...) {
  cat("Model summary (", time_domain_label(x$dt), ")\n", sep = "")
  cat("  states: ", x$dims["states"], ", inputs: ", x$dims["inputs"],
      ", outputs: ", x$dims["outputs"], "\n", sep = "")
  cat("  stability: ", x$stability, "\n", sep = "")
  if (length(x$poles)) {
    cat("  poles: ", paste(format(signif(x$poles, digits)), collapse = ", "),
        "\n", sep = "")
  }
  if (!is.na(x$dims["states"])) {
    cat("  controllability rank: ", x$controllability_rank, " of ",
        x$dims["states"], "\n", sep = "")
    cat("  observability rank: ", x$observability_rank, " of ",
        x$dims["states"], "\n", sep = "")
  }
  if (!is.null(x$dc_gain)) {
    cat("  DC gain:\n")
    print(signif(x$dc_gain, digits))
  }
  invisible(x)
}

#' Select inputs and outputs of a model
#'
#' `sys[i, j]` keeps outputs `i` and inputs `j` (indices or names).
#'
#' @param x A `ctrl_ss` model.
#' @param i Outputs to keep.
#' @param j Inputs to keep.
#' @param ... Unused.
#' @param drop Unused; dimensions are never dropped.
#' @return A `ctrl_ss` model.
#' @examples
#' sys <- ctrl_ss(-diag(2), diag(2), diag(2))
#' sys[1, 2]
#' @export
`[.ctrl_ss` <- function(x, i, j, ..., drop = FALSE) {
  p <- nrow(x$C)
  m <- ncol(x$B)
  oi <- if (missing(i)) seq_len(p) else resolve_index(i, x$output_names, "output")
  ij <- if (missing(j)) seq_len(m) else resolve_index(j, x$input_names, "input")
  ctrl_ss(x$A, x$B[, ij, drop = FALSE], x$C[oi, , drop = FALSE],
          x$D[oi, ij, drop = FALSE], dt = x$dt, state_names = x$state_names,
          input_names = x$input_names[ij], output_names = x$output_names[oi])
}

resolve_index <- function(i, names, what) {
  if (is.character(i)) {
    idx <- match(i, names)
    if (anyNA(idx)) {
      ctrl_abort("Unknown ", what, " name(s): ",
                 paste(i[is.na(idx)], collapse = ", "), ".")
    }
    return(idx)
  }
  if (is.logical(i)) i <- which(rep_len(i, length(names)))
  i <- as.integer(i)
  if (any(i < 1L | i > length(names))) {
    ctrl_abort("Invalid ", what, " index; the model has ", length(names), " ",
               what, "s.")
  }
  i
}

#' Similarity transformation of a state-space model
#'
#' Apply the change of state coordinates \eqn{z = T x}, giving
#' \eqn{(T A T^{-1}, T B, C T^{-1}, D)}.
#'
#' @param sys A `ctrl_ss` model.
#' @param T Nonsingular `n x n` transformation matrix.
#' @return The transformed `ctrl_ss` model (input-output behavior unchanged).
#' @examples
#' sys <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
#' ctrl_transform(sys, matrix(c(1, 1, 0, 1), 2))
#' @export
ctrl_transform <- function(sys, T) {
  sys <- as_ss(sys)
  n <- nrow(sys$A)
  T <- as_real_matrix(T, "T")
  check_dims(T, "T", n, n)
  Ti <- safe_solve(T, diag(n), "transformation matrix `T`")
  ctrl_ss(T %*% sys$A %*% Ti, T %*% sys$B, sys$C %*% Ti, sys$D, dt = sys$dt,
          input_names = sys$input_names, output_names = sys$output_names)
}
