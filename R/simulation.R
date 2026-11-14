#' Simulate a linear model
#'
#' Compute the state and output trajectories of a linear time-invariant
#' model driven by an input signal.
#'
#' **Input format.** `input` is either `NULL` (zero input), a numeric matrix
#' with one row per time point and one column per input (a vector is
#' accepted for single-input models), or a function `function(t)` returning
#' the input vector at time `t` (it is evaluated on the time grid).
#'
#' **Continuous time.** The simulation is exact for the chosen input
#' interpretation. With `hold = "zoh"` the input is held at `u(t_k)` on
#' \eqn{[t_k, t_{k+1})}; with `hold = "foh"` it is linearly interpolated
#' between `u(t_k)` and `u(t_{k+1})`. Each step uses the matrix
#' exponential of the corresponding block matrix (see [ctrl_c2d()]), so
#' there is no integration error. The grid may be irregular; transition
#' matrices are cached for repeated step sizes.
#'
#' **Discrete time.** The recurrence \eqn{x_{k+1} = A x_k + B u_k} is
#' evaluated exactly. `times` must be the sample instants
#' \eqn{t_0 + k\,dt}; `hold` is ignored.
#'
#' Outputs follow \eqn{y_k = C x_k + D u_k} at every time point.
#'
#' @param sys A model (transfer functions are realized with [ctrl_tf2ss()]).
#' @param input Input signal (see details).
#' @param times Strictly increasing time points. For discrete models the
#'   default is `0:steps * dt` with `steps = nrow(input) - 1`.
#' @param x0 Initial state (default zero).
#' @param hold `"zoh"` or `"foh"` input interpretation for continuous models.
#' @return A `ctrl_response` object with `time`, and arrays `x`
#'   (time x states x 1), `y` (time x outputs x 1) and `u` (time x inputs x 1),
#'   plus the signal names. Use [as.data.frame()] for a long data frame and
#'   [plot()] for a plot.
#' @seealso [ctrl_step()], [ctrl_impulse()], [ctrl_initial()]
#' @examples
#' sys <- ctrl_ss(-1, 1, 1)
#' times <- seq(0, 5, by = 0.05)
#' r <- ctrl_simulate(sys, input = sin(times), times = times)
#' head(as.data.frame(r))
#' # Free response of x' = -x from x(0) = 1 equals exp(-t)
#' free <- ctrl_simulate(sys, times = times, x0 = 1)
#' max(abs(free$x[, 1, 1] - exp(-times)))
#' @export
ctrl_simulate <- function(sys, input = NULL, times = NULL, x0 = NULL,
                          hold = c("zoh", "foh")) {
  hold <- match.arg(hold)
  s <- as_ss(sys)
  n <- nrow(s$A)
  m <- ncol(s$B)
  times <- simulation_times(s, times, input)
  nt <- length(times)
  U <- input_matrix(input, times, m)
  x0 <- if (is.null(x0)) rep(0, n) else as.numeric(x0)
  if (length(x0) != n) ctrl_abort("`x0` must have length ", n, ".")
  if (any(!is.finite(x0))) ctrl_abort("`x0` must be finite.")
  X <- matrix(0, nt, n)
  if (n) X[1, ] <- x0
  if (n && nt > 1L) {
    if (s$dt > 0) {
      x <- x0
      for (k in seq_len(nt - 1L)) {
        x <- drop(s$A %*% x + s$B %*% U[k, ])
        X[k + 1L, ] <- x
      }
    } else {
      X <- simulate_continuous(s$A, s$B, U, times, x0, hold)
    }
  }
  Y <- X %*% t(s$C) + U %*% t(s$D)
  new_response(times, X, Y, U, s, type = "simulation")
}

simulation_times <- function(s, times, input) {
  if (is.null(times)) {
    if (s$dt == 0) ctrl_abort("`times` is required for continuous models.")
    nr <- if (is.null(input)) {
      ctrl_abort("Give `times` or an input matrix for discrete simulation.")
    } else if (is.function(input)) {
      ctrl_abort("Give `times` when `input` is a function.")
    } else NROW(input)
    return((seq_len(nr) - 1L) * s$dt)
  }
  if (!is.numeric(times) || length(times) < 1L || any(!is.finite(times))) {
    ctrl_abort("`times` must be a finite numeric vector.")
  }
  if (length(times) > 1L && any(diff(times) <= 0)) {
    ctrl_abort("`times` must be strictly increasing.")
  }
  if (s$dt > 0 && length(times) > 1L) {
    k <- (times - times[1]) / s$dt
    if (max(abs(k - round(k))) > 1e-8 * max(1, max(abs(k))) ||
        any(diff(round(k)) != 1)) {
      ctrl_abort("For a discrete model, `times` must be consecutive sample ",
                 "instants t0 + k * dt with dt = ", s$dt, ".")
    }
  }
  as.numeric(times)
}

input_matrix <- function(input, times, m) {
  nt <- length(times)
  if (is.null(input)) return(matrix(0, nt, m))
  if (is.function(input)) {
    vals <- lapply(times, input)
    U <- do.call(rbind, lapply(vals, as.numeric))
    if (ncol(U) != m) ctrl_abort("`input(t)` must return ", m, " values.")
  } else {
    U <- if (is.null(dim(input))) {
      if (m == 1L) matrix(as.numeric(input), ncol = 1)
      else if (length(input) == m && nt != m) matrix(input, nt, m, byrow = TRUE)
      else ctrl_abort("For a multi-input model, `input` must be a matrix ",
                      "with ", m, " columns.")
    } else as.matrix(input)
    storage.mode(U) <- "double"
    if (nrow(U) == 1L && nt > 1L) U <- matrix(U, nt, m, byrow = TRUE)
  }
  if (nrow(U) != nt) {
    ctrl_abort("`input` must have one row per time point (", nt, "); got ",
               nrow(U), ".")
  }
  if (ncol(U) != m) ctrl_abort("`input` must have ", m, " columns.")
  if (any(!is.finite(U))) ctrl_abort("`input` must be finite.")
  dimnames(U) <- NULL
  U
}

simulate_continuous <- function(A, B, U, times, x0, hold) {
  n <- nrow(A)
  nt <- length(times)
  X <- matrix(0, nt, n)
  X[1, ] <- x0
  h <- diff(times)
  key <- signif(h, 12)
  cache <- new.env(parent = emptyenv())
  x <- x0
  for (k in seq_len(nt - 1L)) {
    kk <- format(key[k], digits = 12)
    mats <- cache[[kk]]
    if (is.null(mats)) {
      mats <- if (hold == "zoh") zoh_matrices(A, B, h[k]) else
        foh_matrices(A, B, h[k])
      cache[[kk]] <- mats
    }
    x <- if (hold == "zoh") {
      drop(mats$A %*% x + mats$B %*% U[k, ])
    } else {
      drop(mats$Phi %*% x + (mats$G1 - mats$G2) %*% U[k, ] +
             mats$G2 %*% U[k + 1L, ])
    }
    X[k + 1L, ] <- x
  }
  X
}

new_response <- function(times, X, Y, U, s, type, cases = NULL) {
  as3 <- function(M) if (length(dim(M)) == 3L) M else
    array(M, c(nrow(M), ncol(M), 1L))
  X <- as3(X)
  Y <- as3(Y)
  U <- as3(U)
  cases <- cases %||% if (dim(Y)[3] == 1L) type else paste0("case", seq_len(dim(Y)[3]))
  structure(
    list(time = times, x = X, y = Y, u = U,
         state_names = s$state_names, input_names = s$input_names,
         output_names = s$output_names, cases = cases, dt = s$dt,
         type = type),
    class = "ctrl_response"
  )
}

#' @export
print.ctrl_response <- function(x, ...) {
  cat("Time response (", x$type, ", ", time_domain_label(x$dt), ")\n", sep = "")
  cat("  ", length(x$time), " time points from ", format(min(x$time)), " to ",
      format(max(x$time)), "\n", sep = "")
  cat("  outputs: ", paste(x$output_names, collapse = ", "), "\n", sep = "")
  if (length(x$cases) > 1L) {
    cat("  cases: ", paste(x$cases, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Long data frame of a time response
#'
#' @param x A `ctrl_response`.
#' @param row.names,optional Unused.
#' @param ... Unused.
#' @return A data frame with columns `time`, `case`, `kind` (`"output"`,
#'   `"state"` or `"input"`), `signal` and `value`.
#' @examples
#' as.data.frame(ctrl_step(ctrl_tf(1, c(1, 1)), times = 0:3))
#' @export
as.data.frame.ctrl_response <- function(x, row.names = NULL, optional = FALSE,
                                        ...) {
  pieces <- list()
  add <- function(arr, names, kind) {
    if (!length(names) || !length(arr)) return()
    for (c in seq_along(x$cases)) for (j in seq_along(names)) {
      pieces[[length(pieces) + 1L]] <<- data.frame(
        time = x$time, case = x$cases[c], kind = kind, signal = names[j],
        value = arr[, j, c], stringsAsFactors = FALSE)
    }
  }
  add(x$y, x$output_names, "output")
  add(x$x, x$state_names, "state")
  add(x$u, x$input_names, "input")
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

# Heuristic final time and step from the pole locations.
auto_time_grid <- function(s, n_points = 1000L) {
  if (s$dt > 0) {
    p <- ctrl_poles(s)
    nz <- p[Mod(p) > 1e-12 & Mod(p) < 1 - 1e-9]
    steps <- if (length(nz)) {
      slow <- max(Mod(nz))
      ceiling(min(5000, max(30, 7 / max(-log(slow), 1e-3))))
    } else 30
    return((0:steps) * s$dt)
  }
  p <- ctrl_poles(s)
  if (!length(p)) return(seq(0, 1, length.out = 101))
  decay <- -Re(p[Re(p) < -1e-9])
  osc <- abs(Im(p[abs(Im(p)) > 1e-9]))
  tf <- if (length(decay)) 7 / min(decay) else 10
  if (any(Re(p) >= -1e-9)) tf <- min(tf, 10 / max(1e-3, max(Mod(p))))
  if (length(osc)) tf <- max(tf, 2 * pi / min(osc) * 1.5)
  fastest <- max(Mod(p), 1e-6)
  h <- min(tf / n_points, 0.1 / fastest)
  npts <- min(20000L, ceiling(tf / h))
  seq(0, tf, length.out = npts + 1L)
}

#' Step, impulse and initial-condition responses
#'
#' * `ctrl_step()`: response to a unit step (times `amplitude`) applied at
#'   time 0 to each selected input in turn, from zero initial state.
#' * `ctrl_impulse()`: impulse response. In continuous time this is the
#'   response to a Dirac impulse, computed as the free response from
#'   \eqn{x(0^+) = B e_j}. The Dirac component \eqn{D\delta(t)} of the
#'   output cannot be represented on a time grid and is omitted (a warning
#'   is given when `D` is nonzero). In discrete time the input is the unit
#'   pulse (Kronecker delta) \eqn{u_0 = 1}, \eqn{u_k = 0} for \eqn{k > 0};
#'   no `1/dt` scaling is applied, so a continuous impulse is *not*
#'   reproduced by a discrete pulse.
#' * `ctrl_initial()`: free response from the initial state `x0`.
#'
#' Default time grids are chosen from the pole locations.
#'
#' @param sys A model.
#' @param times Optional time grid (sample instants for discrete models).
#' @param inputs Input channels (indices or names); default all.
#' @param amplitude Step amplitude.
#' @param x0 Initial state for `ctrl_initial()`.
#' @return A `ctrl_response` with one case per input channel.
#' @examples
#' G <- ctrl_tf(4, c(1, 1.2, 4))
#' st <- ctrl_step(G)
#' ctrl_step_info(st)
#' plot(st)
#' ctrl_impulse(G, times = seq(0, 10, by = 0.05))
#' @export
ctrl_step <- function(sys, times = NULL, inputs = NULL, amplitude = 1) {
  s <- as_ss(sys)
  check_scalar_number(amplitude, "amplitude")
  times <- times %||% auto_time_grid(s)
  times <- simulation_times(s, times, NULL)
  idx <- if (is.null(inputs)) seq_len(ncol(s$B)) else
    resolve_index(inputs, s$input_names, "input")
  responses <- lapply(idx, function(j) {
    U <- matrix(0, length(times), ncol(s$B))
    U[, j] <- amplitude
    ctrl_simulate(s, U, times)
  })
  combine_cases(responses, s, paste("step", s$input_names[idx]), "step")
}

#' @rdname ctrl_step
#' @export
ctrl_impulse <- function(sys, times = NULL, inputs = NULL) {
  s <- as_ss(sys)
  times <- times %||% auto_time_grid(s)
  times <- simulation_times(s, times, NULL)
  idx <- if (is.null(inputs)) seq_len(ncol(s$B)) else
    resolve_index(inputs, s$input_names, "input")
  if (s$dt == 0 && any(s$D[, idx] != 0)) {
    ctrl_warn("The model has direct feedthrough; the Dirac component D * ",
              "delta(t) of the impulse response is not included.")
  }
  responses <- lapply(idx, function(j) {
    if (s$dt > 0) {
      U <- matrix(0, length(times), ncol(s$B))
      U[1, j] <- 1
      ctrl_simulate(s, U, times)
    } else {
      r <- ctrl_simulate(s, NULL, times - times[1], x0 = s$B[, j])
      r$time <- times
      r
    }
  })
  combine_cases(responses, s, paste("impulse", s$input_names[idx]), "impulse")
}

#' @rdname ctrl_step
#' @export
ctrl_initial <- function(sys, x0, times = NULL) {
  s <- as_ss(sys)
  times <- times %||% auto_time_grid(s)
  r <- ctrl_simulate(s, NULL, times, x0 = x0)
  r$type <- "initial"
  r$cases <- "initial"
  r
}

combine_cases <- function(responses, s, cases, type) {
  stack <- function(field) {
    arrs <- lapply(responses, function(r) r[[field]][, , 1, drop = FALSE])
    d <- dim(arrs[[1]])
    out <- array(0, c(d[1], d[2], length(arrs)))
    for (k in seq_along(arrs)) out[, , k] <- arrs[[k]]
    out
  }
  new_response(responses[[1]]$time, stack("x"), stack("y"), stack("u"), s,
               type = type, cases = cases)
}
