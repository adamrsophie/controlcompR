# A static gain (no states, or constant transfer functions) adopts the time
# domain of the other model.
is_static_model <- function(s) {
  if (is_ctrl_ss(s)) return(nrow(s$A) == 0L)
  if (is_ctrl_tf(s)) {
    return(all(vapply(s$den, function(d) sum(d != 0) == 1L && d[length(d)] != 0, TRUE)) &&
             all(vapply(seq_along(s$num), function(k) {
               length(s$num[[k]]) <= length(s$den[[k]]) &&
                 all(s$num[[k]][-length(s$num[[k]])] == 0)
             }, TRUE)))
  }
  FALSE
}

common_dt <- function(s1, s2) {
  if (s1$dt == s2$dt) return(s1$dt)
  if (is_static_model(s1)) return(s2$dt)
  if (is_static_model(s2)) return(s1$dt)
  ctrl_abort("Cannot connect models with different time domains (dt = ",
             s1$dt, " and dt = ", s2$dt, ").")
}

#' Series, parallel and feedback connections
#'
#' * `ctrl_series(sys1, sys2)`: the output of `sys1` drives `sys2`, so the
#'   overall transfer function is \eqn{G_2 G_1} (note the order). The
#'   operator form is `sys2 * sys1`.
#' * `ctrl_parallel(sys1, sys2)`: \eqn{G_1 + G_2} (same inputs, summed
#'   outputs); operator form `sys1 + sys2`.
#' * `ctrl_feedback(sys1, sys2, sign = -1)`: closes the loop
#'   \eqn{u_1 = r + \mathrm{sign} \cdot y_2}, \eqn{u_2 = y_1}, with output
#'   \eqn{y_1}. The default `sign = -1` is negative feedback, giving
#'   \eqn{(I + G_1 G_2)^{-1} G_1}.
#' * `ctrl_append(...)`: block-diagonal (parallel, unconnected) stacking.
#'
#' When both models have direct feedthrough the feedback loop contains an
#' algebraic loop; it is resolved by solving with
#' \eqn{I - \mathrm{sign}\, D_1 D_2}, and an error is raised when this
#' matrix is singular or ill-conditioned (the interconnection is then not
#' well posed).
#'
#' Models must share the time domain; a static gain adopts the time domain
#' of the other model. SISO transfer functions are combined with polynomial
#' arithmetic (without cancellation); all other combinations are carried out
#' in state space.
#'
#' @param sys1,sys2 Models or numeric gains.
#' @param sign Feedback sign, `-1` (negative) or `+1` (positive).
#' @param ... Models to append.
#' @return A model.
#' @examples
#' G <- ctrl_tf(1, c(1, 1))
#' C <- ctrl_tf(c(2, 1), c(1, 0))
#' L <- ctrl_series(C, G)        # G * C
#' T <- ctrl_feedback(L, 1)      # unity negative feedback
#' ctrl_poles(T)
#' @export
ctrl_series <- function(sys1, sys2) {
  sys1 <- as_system(sys1)
  sys2 <- as_system(sys2)
  if (siso_tf_pair(sys1, sys2)) {
    dt <- common_dt(sys1, sys2)
    return(ctrl_tf(poly_mul(sys1$num[[1]], sys2$num[[1]]),
                   poly_mul(sys1$den[[1]], sys2$den[[1]]), dt = dt,
                   input_names = sys1$input_names,
                   output_names = sys2$output_names))
  }
  s1 <- as_ss(sys1)
  s2 <- as_ss(sys2)
  dt <- common_dt(s1, s2)
  if (nrow(s1$C) != ncol(s2$B)) {
    ctrl_abort("Series connection needs ", ncol(s2$B), " outputs from `sys1`; ",
               "it has ", nrow(s1$C), ".")
  }
  n1 <- nrow(s1$A)
  n2 <- nrow(s2$A)
  A <- rbind(cbind(s1$A, matrix(0, n1, n2)), cbind(s2$B %*% s1$C, s2$A))
  B <- rbind(s1$B, s2$B %*% s1$D)
  C <- cbind(s2$D %*% s1$C, s2$C)
  D <- s2$D %*% s1$D
  build_ss(A, B, C, D, dt, s1$input_names, s2$output_names)
}

build_ss <- function(A, B, C, D, dt, input_names, output_names) {
  if (!length(A)) {
    return(ctrl_gain(D, dt = dt, input_names = input_names,
                     output_names = output_names))
  }
  ctrl_ss(A, B, C, D, dt = dt, input_names = input_names,
          output_names = output_names)
}

as_system <- function(x) {
  if (is_ctrl_system(x)) {
    if (inherits(x, "ctrl_zpk")) return(ctrl_zpk2tf(x))
    return(x)
  }
  if (is.numeric(x)) {
    return(if (length(x) == 1L && is.null(dim(x))) ctrl_tf(x, 1) else ctrl_gain(x))
  }
  ctrl_abort("Expected a model or a numeric gain.")
}

siso_tf_pair <- function(a, b) {
  is_ctrl_tf(a) && is_ctrl_tf(b) && length(a$num) == 1L && length(b$num) == 1L
}

#' @rdname ctrl_series
#' @export
ctrl_parallel <- function(sys1, sys2) {
  sys1 <- as_system(sys1)
  sys2 <- as_system(sys2)
  if (siso_tf_pair(sys1, sys2)) {
    dt <- common_dt(sys1, sys2)
    if (identical(sys1$den[[1]], sys2$den[[1]])) {
      return(ctrl_tf(poly_add(sys1$num[[1]], sys2$num[[1]]), sys1$den[[1]],
                     dt = dt, input_names = sys1$input_names,
                     output_names = sys1$output_names))
    }
    return(ctrl_tf(poly_add(poly_mul(sys1$num[[1]], sys2$den[[1]]),
                            poly_mul(sys2$num[[1]], sys1$den[[1]])),
                   poly_mul(sys1$den[[1]], sys2$den[[1]]), dt = dt,
                   input_names = sys1$input_names,
                   output_names = sys1$output_names))
  }
  s1 <- as_ss(sys1)
  s2 <- as_ss(sys2)
  dt <- common_dt(s1, s2)
  if (!identical(dim(s1$D), dim(s2$D))) {
    ctrl_abort("Parallel connection needs models of equal dimensions.")
  }
  build_ss(blkdiag(s1$A, s2$A), rbind(s1$B, s2$B), cbind(s1$C, s2$C),
           s1$D + s2$D, dt, s1$input_names, s1$output_names)
}

#' @rdname ctrl_series
#' @export
ctrl_feedback <- function(sys1, sys2 = 1, sign = -1) {
  if (!sign %in% c(-1, 1)) ctrl_abort("`sign` must be -1 or +1.")
  sys1 <- as_system(sys1)
  sys2 <- as_system(sys2)
  if (siso_tf_pair(sys1, sys2)) {
    dt <- common_dt(sys1, sys2)
    n1 <- sys1$num[[1]]
    d1 <- sys1$den[[1]]
    n2 <- sys2$num[[1]]
    d2 <- sys2$den[[1]]
    # well-posedness: the feedthrough loop gain 1 - sign D1 D2 must not vanish
    # (an improper factor has no finite feedthrough: the check does not
    # apply and the rational result is formed directly)
    dfeed <- function(num, den) {
      num <- strip_leading_zeros(num)
      den <- strip_leading_zeros(den)
      if (length(num) < length(den)) 0 else if (length(num) == length(den))
        num[1] / den[1] else NA_real_
    }
    eloop <- 1 - sign * dfeed(n1, d1) * dfeed(n2, d2)
    if (!is.na(eloop) && abs(eloop) < 1e-12) {
      ctrl_abort("The feedback loop has a singular algebraic loop ",
                 "(1 - sign D1 D2 = 0); the interconnection is not well posed.")
    }
    den <- poly_add(poly_mul(d1, d2), -sign * poly_mul(n1, n2))
    if (all(abs(den) <= 1e-14 * max(abs(c(d1, d2, n1, n2))))) {
      ctrl_abort("The feedback interconnection is not well posed.")
    }
    return(ctrl_tf(poly_mul(n1, d2), den, dt = dt,
                   input_names = sys1$input_names,
                   output_names = sys1$output_names))
  }
  s1 <- as_ss(sys1)
  s2 <- as_ss(sys2)
  dt <- common_dt(s1, s2)
  p1 <- nrow(s1$C)
  m1 <- ncol(s1$B)
  if (ncol(s2$B) != p1 || nrow(s2$C) != m1) {
    ctrl_abort("Feedback needs `sys2` with ", p1, " inputs and ", m1,
               " outputs.")
  }
  n1 <- nrow(s1$A)
  n2 <- nrow(s2$A)
  E <- diag(p1) - sign * s1$D %*% s2$D
  rc <- rcond(E)
  if (rc < 1e-12) {
    ctrl_abort("The feedback loop has an algebraic loop that is singular or ",
               "ill-conditioned (rcond(I - sign D1 D2) = ", signif(rc, 3),
               "); the interconnection is not well posed.")
  }
  rhs <- cbind(s1$C, sign * s1$D %*% s2$C)
  Cy <- if (ncol(rhs)) solve(E, rhs) else matrix(0, p1, 0)
  Dy <- solve(E, s1$D)
  Cu <- sign * cbind(matrix(0, m1, n1), s2$C) + sign * s2$D %*% Cy
  Du <- diag(m1) + sign * s2$D %*% Dy
  A <- blkdiag(s1$A, s2$A) + rbind(s1$B %*% Cu, s2$B %*% Cy)
  B <- rbind(s1$B %*% Du, s2$B %*% Dy)
  build_ss(A, B, Cy, Dy, dt, s1$input_names, s1$output_names)
}

#' @rdname ctrl_series
#' @export
ctrl_append <- function(...) {
  ss <- lapply(list(...), function(x) as_ss(as_system(x)))
  # static gains adopt the time domain of the dynamic blocks
  dyn <- unique(vapply(Filter(function(x) nrow(x$A) > 0L, ss), `[[`, 0, "dt"))
  if (length(dyn) > 1L) ctrl_abort("Cannot append models with different time domains.")
  dt <- if (length(dyn)) dyn else ss[[1]]$dt
  A <- blkdiag(lapply(ss, `[[`, "A"))
  B <- blkdiag(lapply(ss, `[[`, "B"))
  C <- blkdiag(lapply(ss, `[[`, "C"))
  D <- blkdiag(lapply(ss, `[[`, "D"))
  inn <- unlist(lapply(ss, `[[`, "input_names"))
  out <- unlist(lapply(ss, `[[`, "output_names"))
  if (anyDuplicated(inn)) inn <- paste0("u", seq_along(inn))
  if (anyDuplicated(out)) out <- paste0("y", seq_along(out))
  build_ss(A, B, C, D, dt, inn, out)
}

#' Arithmetic on models
#'
#' `+` is a parallel connection, `-` subtracts (or negates), `*` is a series
#' connection with MATLAB semantics (`G2 * G1` applies `G1` first), `/`
#' divides by a scalar or (SISO transfer functions) multiplies by the
#' inverse, and `^` raises a model to a nonnegative integer power (negative
#' powers for SISO transfer functions).
#'
#' @param e1,e2 Models or numbers.
#' @return A model.
#' @examples
#' s <- ctrl_tf("s")
#' G <- 1 / (s + 1)
#' 2 * G + G^2
#' @name ctrl_arith
#' @export
Ops.ctrl_system <- function(e1, e2) {
  op <- .Generic
  if (missing(e2)) {
    if (op == "-") return(negate_system(as_system(e1)))
    if (op == "+") return(e1)
    ctrl_abort("Unsupported unary operator `", op, "`.")
  }
  switch(op,
    "+" = ctrl_parallel(e1, e2),
    "-" = ctrl_parallel(e1, negate_system(as_system(e2))),
    "*" = {
      # a numeric scalar scales every channel; matrices and models connect
      # in series
      if (is_scalar_number(e1)) scale_system(as_system(e2), e1) else
        if (is_scalar_number(e2)) scale_system(as_system(e1), e2) else
          ctrl_series(e2, e1)
    },
    "/" = {
      if (is_scalar_number(e2)) {
        if (e2 == 0) ctrl_abort("Division of a model by zero.")
        scale_system(as_system(e1), 1 / e2)
      } else if (is.numeric(e2)) {
        ctrl_series(as_system(e1), as_system(1 / e2))
      } else {
        ctrl_series(ctrl_inv(e2), e1)
      }
    },
    "^" = {
      if (!is.numeric(e2) || length(e2) != 1L || e2 != round(e2)) {
        ctrl_abort("Model powers must be integers.")
      }
      base <- as_system(e1)
      if (e2 < 0) {
        base <- ctrl_inv(base)
        e2 <- -e2
      }
      out <- if (is_ctrl_tf(base)) ctrl_tf(1, 1, dt = base$dt) else
        ctrl_gain(diag(nrow(as_ss(base)$C)), dt = base$dt)
      for (i in seq_len(e2)) out <- ctrl_series(out, base)
      out
    },
    "==" = ctrl_abort("Use all.equal() or compare frequency responses to ",
                      "compare models."),
    ctrl_abort("Unsupported operator `", op, "` for models.")
  )
}

is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.null(dim(x)) && is.finite(x)
}

negate_system <- function(s) scale_system(s, -1)

# Multiply every channel of a model by the scalar k (names and time domain
# are kept).
scale_system <- function(s, k) {
  if (is_ctrl_tf(s)) {
    s$num <- matrix(lapply(s$num, function(p) k * p), nrow(s$num), ncol(s$num))
    return(s)
  }
  s <- as_ss(s)
  s$C <- k * s$C
  s$D <- k * s$D
  s
}

#' Inverse of a model
#'
#' SISO transfer functions are inverted by swapping numerator and
#' denominator (the result may be improper). State-space models need a
#' square, invertible feedthrough `D`; the inverse is
#' \eqn{(A - B D^{-1} C, B D^{-1}, -D^{-1} C, D^{-1})}.
#'
#' @param sys A model.
#' @return The inverse model.
#' @examples
#' ctrl_inv(ctrl_tf(c(1, 2), c(1, 3)))
#' @export
ctrl_inv <- function(sys) {
  sys <- as_system(sys)
  if (is_ctrl_tf(sys) && length(sys$num) == 1L) {
    if (all(sys$num[[1]] == 0)) ctrl_abort("Cannot invert a zero model.")
    return(ctrl_tf(sys$den[[1]], sys$num[[1]], dt = sys$dt,
                   input_names = sys$output_names,
                   output_names = sys$input_names))
  }
  s <- as_ss(sys)
  if (nrow(s$D) != ncol(s$D)) ctrl_abort("Only square models can be inverted.")
  Di <- safe_solve(s$D, diag(nrow(s$D)), "feedthrough matrix D")
  build_ss(s$A - s$B %*% Di %*% s$C, s$B %*% Di, -Di %*% s$C, Di, s$dt,
           s$output_names, s$input_names)
}

#' Summing junction
#'
#' Create a static gain block from an equation such as `"e = r - y"`. Use it
#' with [ctrl_connect()] to build block diagrams by signal name.
#'
#' @param expr A string `"out = +a - b + c"`.
#' @param dt Time domain (a static block adopts the domain of the others).
#' @return A `ctrl_ss` static gain with named inputs and output.
#' @examples
#' ctrl_sumblock("e = r - y")
#' @export
ctrl_sumblock <- function(expr, dt = NULL) {
  if (!is.character(expr) || length(expr) != 1L || !grepl("=", expr)) {
    ctrl_abort("`expr` must be a string like \"e = r - y\".")
  }
  parts <- strsplit(expr, "=", fixed = TRUE)[[1]]
  lhs <- trimws(parts[1])
  rhs <- gsub("\\s+", "", parts[2])
  if (!grepl("^[+-]", rhs)) rhs <- paste0("+", rhs)
  toks <- regmatches(rhs, gregexpr("[+-][^+-]+", rhs))[[1]]
  sgn <- ifelse(substr(toks, 1, 1) == "-", -1, 1)
  nm <- substring(toks, 2)
  ctrl_gain(matrix(sgn, 1), dt = dt, input_names = nm, output_names = lhs)
}

#' Connect models by signal names
#'
#' Build an interconnection of several models whose input and output names
#' define the wiring, as in a block diagram. Every block input whose name
#' equals a block output name is driven by that output; the remaining block
#' inputs must be listed in `inputs` and become external inputs. Use
#' [ctrl_sumblock()] for sums of signals. Output names must be unique across
#' blocks.
#'
#' Algebraic loops through feedthrough terms are resolved exactly; an
#' error is raised when the loop matrix is singular (ill-posed diagram).
#'
#' @param ... Models with named signals.
#' @param inputs Names of the external inputs.
#' @param outputs Names of the block outputs to expose.
#' @return A `ctrl_ss` model from `inputs` to `outputs`.
#' @examples
#' G <- ctrl_tf(1, c(1, 1), input_names = "u", output_names = "y")
#' K <- ctrl_tf(5, 1, input_names = "e", output_names = "u")
#' S <- ctrl_sumblock("e = r - y")
#' T <- ctrl_connect(G, K, S, inputs = "r", outputs = "y")
#' ctrl_dcgain(T)  # 5 / 6
#' @export
ctrl_connect <- function(..., inputs, outputs) {
  blocks <- lapply(list(...), function(b) as_ss(as_system(b)))
  dts <- unique(vapply(Filter(function(b) nrow(b$A) > 0L, blocks),
                       `[[`, 0, "dt"))
  if (length(dts) > 1L) ctrl_abort("All blocks must share the time domain.")
  dt <- if (length(dts)) dts else 0
  A <- blkdiag(lapply(blocks, `[[`, "A"))
  B <- blkdiag(lapply(blocks, `[[`, "B"))
  C <- blkdiag(lapply(blocks, `[[`, "C"))
  D <- blkdiag(lapply(blocks, `[[`, "D"))
  in_names <- unlist(lapply(blocks, `[[`, "input_names"))
  out_names <- unlist(lapply(blocks, `[[`, "output_names"))
  if (anyDuplicated(out_names)) {
    ctrl_abort("Block output names must be unique; duplicated: ",
               paste(unique(out_names[duplicated(out_names)]), collapse = ", "))
  }
  inputs <- as.character(inputs)
  outputs <- as.character(outputs)
  F <- matrix(0, length(in_names), length(out_names))
  G <- matrix(0, length(in_names), length(inputs))
  for (i in seq_along(in_names)) {
    j <- match(in_names[i], out_names)
    if (!is.na(j)) {
      F[i, j] <- 1
    } else {
      k <- match(in_names[i], inputs)
      if (is.na(k)) {
        ctrl_abort("Block input \"", in_names[i], "\" is neither a block ",
                   "output nor listed in `inputs`.")
      }
      G[i, k] <- 1
    }
  }
  sel <- match(outputs, out_names)
  if (anyNA(sel)) {
    ctrl_abort("Unknown output name(s): ",
               paste(outputs[is.na(sel)], collapse = ", "), ".")
  }
  E <- diag(length(out_names)) - D %*% F
  rc <- rcond(E)
  if (rc < 1e-12) {
    ctrl_abort("The block diagram has a singular algebraic loop (rcond ",
               signif(rc, 3), ").")
  }
  Cy <- solve(E, C)
  Dy <- solve(E, D %*% G)
  Anew <- A + B %*% F %*% Cy
  Bnew <- B %*% (F %*% Dy + G)
  build_ss(Anew, Bnew, Cy[sel, , drop = FALSE], Dy[sel, , drop = FALSE], dt,
           inputs, outputs)
}
