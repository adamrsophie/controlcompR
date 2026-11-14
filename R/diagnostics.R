#' Step-response characteristics
#'
#' Compute the usual transient characteristics for each output and input
#' case of a step response. Definitions (with \eqn{y_0} the value before
#' the step and \eqn{y_\infty} the final value; when the model is known,
#' a direct feedthrough jump at \eqn{t = 0} is not counted in \eqn{y_0}):
#'
#' * rise time: time to go from `rise[1]` to `rise[2]` (default 10 to 90
#'   percent) of \eqn{y_\infty - y_0}, with linear interpolation;
#' * settling time: the last time the response is outside the band
#'   \eqn{|y - y_\infty| \le} `settle` \eqn{\cdot |y_\infty - y_0|}
#'   (default 2 percent); `NA` if it never enters the band for good on the
#'   simulated horizon;
#' * overshoot and undershoot: percent of \eqn{|y_\infty - y_0|};
#' * peak and peak time: largest \eqn{|y - y_0|}.
#'
#' When `sys` (the model) is supplied the final value is its DC gain times
#' the step amplitude (exact, for stable models). Otherwise the last sample
#' is used and `settled` reports whether the response varies by less than
#' the band over the last 5 percent of the horizon. Characteristics depend
#' on these conventions, so state them when reporting numbers.
#'
#' @param response A `ctrl_response` (typically from [ctrl_step()]) or a
#'   model, whose step response is then computed.
#' @param sys Optional model for the exact final value.
#' @param settle Relative settling band.
#' @param rise Rise-time fractions.
#' @return A data frame with one row per output and case.
#' @examples
#' G <- ctrl_tf(4, c(1, 1.2, 4))
#' ctrl_step_info(G)
#' @export
ctrl_step_info <- function(response, sys = NULL, settle = 0.02,
                           rise = c(0.1, 0.9)) {
  if (is_ctrl_system(response)) {
    sys <- sys %||% response
    response <- ctrl_step(response)
  }
  if (!inherits(response, "ctrl_response")) {
    ctrl_abort("`response` must be a `ctrl_response` or a model.")
  }
  t <- response$time
  dc <- NULL
  amp <- 1
  if (!is.null(sys)) {
    dc <- tryCatch(ctrl_dcgain(sys), error = function(e) NULL)
    if (!ctrl_is_stable(sys)) dc <- NULL
  }
  Dm <- if (!is.null(sys)) tryCatch(as_ss(as_system(sys))$D,
                                     error = function(e) NULL) else NULL
  rows <- list()
  for (c in seq_along(response$cases)) for (i in seq_along(response$output_names)) {
    y <- response$y[, i, c]
    y0 <- y[1]
    # a model with direct feedthrough jumps at t = 0; the value before the
    # step is y[1] - D u[1]
    if (!is.null(Dm)) y0 <- y0 - sum(Dm[i, ] * response$u[1, , c])
    j <- which(response$u[length(t), , c] != 0)[1]
    if (!is.null(dc) && !is.na(j)) {
      yinf <- y0 + dc[i, j] * response$u[length(t), j, c]
      settled <- TRUE
    } else {
      yinf <- y[length(y)]
      tail_idx <- t >= t[1] + 0.95 * (t[length(t)] - t[1])
      settled <- diff(range(y[tail_idx])) <= settle * abs(yinf - y0)
    }
    dy <- yinf - y0
    rows[[length(rows) + 1L]] <- step_metrics(t, y, y0, yinf, dy, settle, rise,
                                              response$output_names[i],
                                              response$cases[c], settled)
  }
  do.call(rbind, rows)
}

step_metrics <- function(t, y, y0, yinf, dy, settle, rise, out, case, settled) {
  na <- NA_real_
  if (dy == 0) {
    return(data.frame(output = out, case = case, rise_time = na,
                      settling_time = na, overshoot = na, undershoot = na,
                      peak = max(abs(y - y0)), peak_time = t[which.max(abs(y - y0))],
                      final_value = yinf, settled = settled))
  }
  frac <- (y - y0) / dy
  cross <- function(level) {
    k <- which(frac >= level)[1]
    if (is.na(k)) return(na)
    if (k == 1L) return(t[1])
    t[k - 1L] + (level - frac[k - 1L]) / (frac[k] - frac[k - 1L]) * (t[k] - t[k - 1L])
  }
  tr <- cross(rise[2]) - cross(rise[1])
  out_band <- which(abs(y - yinf) > settle * abs(dy))
  ts <- if (!length(out_band)) t[1] else if (max(out_band) == length(t)) na else
    t[max(out_band) + 1L]
  os <- max(0, max(frac) - 1) * 100
  us <- max(0, -min(frac)) * 100
  k <- which.max(abs(y - y0))
  data.frame(output = out, case = case, rise_time = tr, settling_time = ts,
             overshoot = os, undershoot = us, peak = abs(y[k] - y0) + 0 * y0,
             peak_time = t[k], final_value = yinf, settled = settled)
}

#' Tracking and control performance metrics
#'
#' Integral criteria for a tracking or regulation experiment on a time
#' grid (trapezoidal rule in continuous time, sums times `dt` for sampled
#' data):
#' IAE \eqn{\int |e|}, ISE \eqn{\int e^2}, ITAE \eqn{\int t|e|}, RMSE
#' \eqn{\sqrt{\mathrm{mean}(e^2)}}, control effort \eqn{\int u^2}, maximum
#' \eqn{|u|} and, when limits are given, the largest constraint violation.
#'
#' @param time Time vector, with one point per row of `y` (one extra final
#'   point, as returned by simulations that include the final state, is
#'   dropped).
#' @param y Output matrix or vector.
#' @param r Reference (same shape as `y`, or a scalar).
#' @param u Optional control signal.
#' @param u_min,u_max Optional input limits for violation reporting, a
#'   scalar or one value per input (column of `u`).
#' @param discrete Use rectangle sums instead of the trapezoidal rule.
#' @return A one-row data frame.
#' @examples
#' t <- seq(0, 10, by = 0.01)
#' ctrl_performance(t, 1 - exp(-t), 1)
#' @export
ctrl_performance <- function(time, y, r = 0, u = NULL, u_min = -Inf,
                             u_max = Inf, discrete = FALSE) {
  y <- as.matrix(y)
  time <- as.numeric(time)
  # simulations often return one more time point (the final state) than
  # output samples
  if (length(time) == nrow(y) + 1L) time <- time[seq_len(nrow(y))]
  if (length(time) != nrow(y)) {
    ctrl_abort("`time` has ", length(time), " points but `y` has ", nrow(y),
               " rows.")
  }
  r <- if (length(r) == 1L) matrix(r, nrow(y), ncol(y)) else as.matrix(r)
  if (nrow(r) == nrow(y) + 1L) r <- r[seq_len(nrow(y)), , drop = FALSE]
  if (!identical(dim(r), dim(y))) ctrl_abort("`r` must match the shape of `y`.")
  e <- r - y
  integ <- function(f) {
    f <- as.matrix(f)
    if (discrete || length(time) < 2L) {
      h <- if (length(time) > 1L) stats::median(diff(time)) else 1
      sum(f) * h
    } else {
      dtv <- diff(time)
      sum(colSums((f[-1, , drop = FALSE] + f[-nrow(f), , drop = FALSE]) / 2 * dtv))
    }
  }
  out <- data.frame(IAE = integ(abs(e)), ISE = integ(e^2),
                    ITAE = integ(abs(e) * time), RMSE = sqrt(mean(e^2)))
  if (!is.null(u)) {
    u <- as.matrix(u)
    if (nrow(u) != nrow(y)) {
      ctrl_abort("`u` has ", nrow(u), " rows but `y` has ", nrow(y), ".")
    }
    lim <- function(v, name) {
      if (!length(v) %in% c(1L, ncol(u))) {
        ctrl_abort("`", name, "` must have length 1 or ", ncol(u), ".")
      }
      matrix(rep_len(v, ncol(u)), nrow(u), ncol(u), byrow = TRUE)
    }
    out$control_effort <- integ(u^2)
    out$max_abs_u <- max(abs(u))
    out$max_violation <- max(0, max(u - lim(u_max, "u_max")),
                             max(lim(u_min, "u_min") - u))
  }
  out
}

#' Compare controllers on a common scenario
#'
#' Run several closed-loop simulation functions on the same plant and
#' scenario and collect their performance metrics in one table. Each
#' element of `controllers` is a function of no arguments that returns a
#' list (or data frame) with `time`, `y`, and optionally `u` and `r`; this
#' keeps the comparison independent of the controller type (PID, LQR, MPC
#' or user code).
#'
#' @param controllers Named list of simulation functions.
#' @param reference Reference used when a result has no `r`.
#' @param u_min,u_max Input limits for violation reporting.
#' @return A data frame with one row per controller.
#' @examples
#' plant <- ctrl_tf(1, c(1, 1, 0))
#' sims <- list(
#'   pid_slow = function() ctrl_simulate_pid(plant, ctrl_pid(1, 0.1, 1,
#'     derivative_filter = 10, dt = 0.05), steps = 300),
#'   pid_fast = function() ctrl_simulate_pid(plant, ctrl_pid(4, 0.5, 2,
#'     derivative_filter = 10, dt = 0.05), steps = 300)
#' )
#' ctrl_compare_controllers(sims, reference = 1)
#' @export
ctrl_compare_controllers <- function(controllers, reference = 0, u_min = -Inf,
                                     u_max = Inf) {
  if (!is.list(controllers) || is.null(names(controllers))) {
    ctrl_abort("`controllers` must be a named list of functions.")
  }
  rows <- lapply(names(controllers), function(nm) {
    res <- controllers[[nm]]()
    r <- res$r %||% reference
    perf <- ctrl_performance(res$time, res$y, r, u = res$u, u_min = u_min,
                             u_max = u_max)
    cbind(controller = nm, perf, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
