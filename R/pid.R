#' PID controllers
#'
#' Define a PID controller with the features needed in practice:
#' derivative filtering, setpoint weighting, output saturation, anti-windup
#' and bumpless transfer. The control law (parallel form) is
#' \deqn{u = k_p (b\,r - y) + k_i \int (r - y)\,dt + k_d \frac{d}{dt}
#' (c\,r - y)_f,}
#' where the derivative acts on the filtered signal
#' \eqn{(c r - y)_f = (c r - y) / (T_f s + 1)} and `b`, `c` are setpoint
#' weights (`c = 0` avoids derivative kicks on setpoint steps). With
#' `form = "standard"` the gains are given as \eqn{k_p}, \eqn{T_i} and
#' \eqn{T_d} with \eqn{k_i = k_p / T_i} and \eqn{k_d = k_p T_d}.
#'
#' The derivative filter time constant `Tf` is required whenever
#' \eqn{k_d \ne 0}: an unfiltered derivative is improper and amplifies
#' measurement noise without bound. A common choice is
#' \eqn{T_f = T_d / N} with \eqn{N} between 2 and 20 (`derivative_filter =
#' N`).
#'
#' **Discrete implementation** (used by [ctrl_simulate_pid()] and
#' [ctrl_pid_update()]) with sample interval `dt`: backward Euler for the
#' derivative filter, forward Euler for the integral by default
#' (`integral_method = "forward"`; `"backward"` and `"tustin"` are also
#' available), and saturation to `[u_min, u_max]`. Anti-windup options:
#' `"backcalc"` (back-calculation, the integrator is driven by
#' \eqn{k_t (u_{sat} - u)} with tracking gain `kt`, default
#' \eqn{1/\sqrt{T_i T_d}} or \eqn{1/T_i}), `"clamp"` (conditional
#' integration: the integrator is frozen while saturated in the direction
#' of the error) or `"none"`. Bumpless transfer from manual to automatic
#' mode initializes the integrator so that the first automatic output
#' equals the last manual output.
#'
#' @param kp Proportional gain.
#' @param ki Integral gain (parallel form) or integral time `Ti`
#'   (standard form; `Inf` disables integral action).
#' @param kd Derivative gain (parallel form) or derivative time `Td`
#'   (standard form).
#' @param Tf Derivative filter time constant.
#' @param derivative_filter Alternative to `Tf`: the ratio `N` with
#'   `Tf = Td / N`.
#' @param b,c Setpoint weights for the proportional and derivative terms.
#' @param u_min,u_max Output saturation limits.
#' @param antiwindup `"backcalc"`, `"clamp"` or `"none"`.
#' @param kt Back-calculation tracking gain.
#' @param form `"parallel"` or `"standard"`.
#' @param dt Sample interval of the discrete implementation (`NULL` for a
#'   continuous description; required for simulation).
#' @param integral_method Discretization of the integral term.
#' @return An object of class `c("ctrl_pid", "ctrl_controller")`.
#' @seealso [ctrl_pid_tf()], [ctrl_simulate_pid()], [ctrl_pid_tune()]
#' @examples
#' pid <- ctrl_pid(kp = 2, ki = 1, kd = 0.5, derivative_filter = 10,
#'                 u_min = -1, u_max = 1, dt = 0.01)
#' pid
#' ctrl_pid_tf(pid)
#' @export
ctrl_pid <- function(kp = 1, ki = 0, kd = 0, Tf = NULL, derivative_filter = NULL,
                     b = 1, c = 0, u_min = -Inf, u_max = Inf,
                     antiwindup = c("backcalc", "clamp", "none"), kt = NULL,
                     form = c("parallel", "standard"), dt = NULL,
                     integral_method = c("forward", "backward", "tustin")) {
  antiwindup <- match.arg(antiwindup)
  form <- match.arg(form)
  integral_method <- match.arg(integral_method)
  for (nm in c("kp", "b", "c")) check_scalar_number(get(nm), nm)
  if (!is.numeric(ki) || length(ki) != 1L || is.na(ki)) ctrl_abort("`ki` must be a number.")
  check_scalar_number(kd, "kd")
  if (form == "standard") {
    Ti <- ki
    Td <- kd
    if (Ti <= 0) ctrl_abort("The integral time must be positive (Inf disables it).")
    ki <- if (is.infinite(Ti)) 0 else kp / Ti
    kd <- kp * Td
  } else if (!is.finite(ki)) {
    ctrl_abort("`ki` must be finite in parallel form.")
  }
  if (!is.numeric(u_min) || !is.numeric(u_max) || u_min >= u_max) {
    ctrl_abort("`u_min` must be below `u_max`.")
  }
  if (kd != 0) {
    if (is.null(Tf)) {
      if (is.null(derivative_filter)) {
        ctrl_abort("A derivative term needs a filter: give `Tf` or ",
                   "`derivative_filter` (the ratio N with Tf = Td / N).")
      }
      check_scalar_number(derivative_filter, "derivative_filter", positive = TRUE)
      Td <- kd / kp
      if (kp == 0) ctrl_abort("`derivative_filter` needs kp != 0; give `Tf`.")
      Tf <- abs(Td) / derivative_filter
    }
    check_scalar_number(Tf, "Tf", positive = TRUE)
  } else {
    Tf <- Tf %||% 0
  }
  if (is.null(kt)) {
    kt <- if (ki != 0 && kp != 0) {
      Ti <- kp / ki
      if (kd != 0) 1 / sqrt(abs(Ti * kd / kp)) else 1 / abs(Ti)
    } else 1
  }
  dt <- if (is.null(dt)) 0 else match_time_domain(dt)
  if (dt > 0 && kd != 0 && Tf < dt / 10) {
    ctrl_warn("The derivative filter time constant is much shorter than the ",
              "sample interval; the derivative term will amplify noise.")
  }
  structure(list(kp = kp, ki = ki, kd = kd, Tf = Tf, b = b, c = c,
                 u_min = u_min, u_max = u_max, antiwindup = antiwindup, kt = kt,
                 dt = dt, integral_method = integral_method),
            class = c("ctrl_pid", "ctrl_controller"))
}

#' @export
print.ctrl_pid <- function(x, digits = 4, ...) {
  cat("PID controller (", if (x$dt > 0) paste0("discrete, dt = ", x$dt) else
    "continuous description", ")\n", sep = "")
  cat("  kp = ", signif(x$kp, digits), ", ki = ", signif(x$ki, digits),
      ", kd = ", signif(x$kd, digits), ", Tf = ", signif(x$Tf, digits), "\n",
      sep = "")
  cat("  setpoint weights b = ", x$b, ", c = ", x$c, "\n", sep = "")
  cat("  saturation [", x$u_min, ", ", x$u_max, "], anti-windup: ",
      x$antiwindup, "\n", sep = "")
  invisible(x)
}

#' Linear transfer functions of a PID controller
#'
#' Return the linear part of the controller (ignoring saturation) as the
#' transfer function from the error \eqn{e = r - y} to `u` (this one-input
#' form corresponds to setpoint weights `b = c = 1`),
#' \eqn{C(s) = k_p + k_i/s + k_d s/(T_f s + 1)}, or, with
#' `two_degree = TRUE`, as a two-input model from `(r, y)` to `u` that
#' includes the setpoint weights. When the controller has a sample
#' interval, the discrete equivalent of the implemented difference
#' equations is returned.
#'
#' @param pid A `ctrl_pid` object.
#' @param two_degree Return the two-input form.
#' @return A `ctrl_tf` model.
#' @examples
#' ctrl_pid_tf(ctrl_pid(kp = 1, ki = 2))
#' @export
ctrl_pid_tf <- function(pid, two_degree = FALSE) {
  if (!inherits(pid, "ctrl_pid")) ctrl_abort("`pid` must be a `ctrl_pid` object.")
  dt <- pid$dt
  if (dt == 0) {
    I <- ctrl_tf(pid$ki, c(1, 0))
    D <- if (pid$kd != 0) ctrl_tf(c(pid$kd, 0), c(pid$Tf, 1)) else ctrl_tf(0, 1)
  } else {
    I <- switch(pid$integral_method,
      forward = ctrl_tf(pid$ki * dt, c(1, -1), dt = dt),
      backward = ctrl_tf(c(pid$ki * dt, 0), c(1, -1), dt = dt),
      tustin = ctrl_tf(c(pid$ki * dt / 2, pid$ki * dt / 2), c(1, -1), dt = dt))
    D <- if (pid$kd != 0) {
      # backward Euler: d_k = Tf/(Tf+dt) d_{k-1} + kd/(Tf+dt) (e_k - e_{k-1})
      a <- pid$Tf / (pid$Tf + dt)
      g <- pid$kd / (pid$Tf + dt)
      ctrl_tf(c(g, -g), c(1, -a), dt = dt)
    } else ctrl_tf(0, 1, dt = dt)
  }
  if (!two_degree) {
    P <- ctrl_tf(pid$kp, 1, dt = if (dt > 0) dt else NULL)
    out <- P + I + D
    out$input_names <- "e"
    out$output_names <- "u"
    return(out)
  }
  Pr <- pid$kp * pid$b
  num <- matrix(list(), 1, 2)
  den <- matrix(list(), 1, 2)
  r_tf <- ctrl_tf(Pr, 1, dt = if (dt > 0) dt else NULL) + I + pid$c * D
  y_tf <- -(ctrl_tf(pid$kp, 1, dt = if (dt > 0) dt else NULL) + I + D)
  num[[1, 1]] <- r_tf$num[[1]]
  den[[1, 1]] <- r_tf$den[[1]]
  num[[1, 2]] <- y_tf$num[[1]]
  den[[1, 2]] <- y_tf$den[[1]]
  ctrl_tf(num, den, dt = if (dt > 0) dt else NULL, input_names = c("r", "y"),
          output_names = "u")
}

#' One step of a discrete PID controller
#'
#' Advance the discrete PID implementation by one sample. The state holds
#' the integrator, the filtered derivative, the previous derivative input
#' and the mode.
#'
#' @param pid A `ctrl_pid` object with a positive `dt`.
#' @param state Controller state from the previous call (`NULL` to
#'   initialize).
#' @param r,y Reference and measurement at this sample.
#' @param manual Optional manual output; when given, the controller is in
#'   manual mode and tracks this value for bumpless transfer.
#' @return A list with `u` (saturated output), `u_unsat`, the parts `P`,
#'   `I`, `D`, `saturated` and the updated `state`.
#' @examples
#' pid <- ctrl_pid(kp = 1, ki = 1, dt = 0.1)
#' st <- NULL
#' for (k in 1:3) {
#'   out <- ctrl_pid_update(pid, st, r = 1, y = 0)
#'   st <- out$state
#' }
#' out$u
#' @export
ctrl_pid_update <- function(pid, state = NULL, r, y, manual = NULL) {
  dt <- pid$dt
  if (dt <= 0) ctrl_abort("The PID controller needs a sample interval `dt`.")
  if (is.null(state)) {
    state <- list(I = 0, D = 0, ed_prev = pid$c * r - y, e_prev = r - y,
                  manual_prev = FALSE)
  }
  e <- r - y
  P <- pid$kp * (pid$b * r - y)
  ed <- pid$c * r - y
  D <- if (pid$kd != 0) {
    a <- pid$Tf / (pid$Tf + dt)
    a * state$D + pid$kd / (pid$Tf + dt) * (ed - state$ed_prev)
  } else 0
  I <- state$I
  if (pid$integral_method == "backward") I <- I + pid$ki * dt * e
  if (pid$integral_method == "tustin") I <- I + pid$ki * dt * (e + state$e_prev) / 2
  if (!is.null(manual)) {
    # track the manual output: the integrator absorbs the difference so the
    # first automatic output equals the manual one when r and y are unchanged
    # (it anticipates that sample's integral increment and derivative decay)
    u <- min(max(manual, pid$u_min), pid$u_max)
    inc <- if (pid$integral_method == "forward") 0 else pid$ki * dt * e
    D_next <- if (pid$kd != 0) pid$Tf / (pid$Tf + dt) * D else 0
    I_next <- u - P - D_next - inc
    return(list(u = u, u_unsat = u, P = P, I = I, D = D, saturated = FALSE,
                state = list(I = I_next, D = D, ed_prev = ed, e_prev = e,
                             manual_prev = TRUE)))
  }
  v <- P + I + D
  u <- min(max(v, pid$u_min), pid$u_max)
  sat <- u != v
  I_next <- I
  grow <- switch(pid$integral_method,
                 forward = pid$ki * dt * e,
                 backward = 0,
                 tustin = 0)
  if (pid$antiwindup == "backcalc") {
    I_next <- I + grow + pid$kt * dt * (u - v)
  } else if (pid$antiwindup == "clamp") {
    pushing <- (v > pid$u_max && e * pid$ki > 0) || (v < pid$u_min && e * pid$ki < 0)
    if (!pushing) I_next <- I + grow
    if (pushing && pid$integral_method != "forward") I_next <- state$I
  } else {
    I_next <- I + grow
  }
  list(u = u, u_unsat = v, P = P, I = I, D = D, saturated = sat,
       state = list(I = I_next, D = D, ed_prev = ed, e_prev = e,
                    manual_prev = FALSE))
}

# Solve u = sat(alpha - beta * (y0 + Dp * u)) for one sample. The PID
# update is affine in the measurement before saturation, so alpha and beta
# come from two evaluations.
pid_loop_solve <- function(pid, st, r, y0, Dp) {
  v0 <- ctrl_pid_update(pid, st, r, y0)$u_unsat
  v1 <- ctrl_pid_update(pid, st, r, y0 + 1)$u_unsat
  beta <- v0 - v1
  kk <- beta * Dp
  if (!(1 + kk > 1e-12)) {
    ctrl_abort("The PID loop through the plant feedthrough is not well posed ",
               "(1 + beta Dp = ", signif(1 + kk, 4), ").")
  }
  # unsaturated solution of u = v0 - kk u, else the active bound
  u <- v0 / (1 + kk)
  min(max(u, pid$u_min), pid$u_max)
}

#' Closed-loop simulation with a discrete PID controller
#'
#' Simulate a plant under the discrete PID implementation, including
#' saturation, anti-windup and optional manual mode. A continuous plant is
#' discretized exactly with a zero-order hold at the controller sample
#' interval (the controller output is held between samples, as with a
#' digital controller and DAC); a discrete plant must have the controller's
#' sample interval. The measured output is the first plant output unless
#' `output` is given.
#'
#' With direct feedthrough \eqn{D_p} from the plant input to the measured
#' output, the measurement depends on the input of the same sample. The
#' unsaturated controller output is affine in the measurement,
#' \eqn{v = \alpha - \beta y}, so the loop
#' \eqn{u = \mathrm{sat}(\alpha - \beta (y_0 + D_p u))} is solved exactly;
#' it has a unique solution when \eqn{1 + \beta D_p > 0} and an error is
#' raised otherwise.
#'
#' @param plant SISO plant (or a plant whose input `input` is driven).
#' @param pid A `ctrl_pid` object with positive `dt`.
#' @param reference Reference signal: a number, a vector of length `steps`
#'   or a function of time.
#' @param steps Number of samples.
#' @param disturbance Optional input disturbance (added to the plant input),
#'   same formats as `reference`.
#' @param noise Optional measurement noise, same formats.
#' @param manual Optional vector of manual outputs (`NA` where automatic).
#' @param x0 Initial plant state.
#' @param output Index of the measured plant output.
#' @return A data frame with `time`, `r`, `y` (measured), `y_true`, `u`,
#'   `u_unsat`, `P`, `I`, `D`, `saturated` and `mode`, with class
#'   `ctrl_pid_simulation`.
#' @examples
#' plant <- ctrl_tf(1, c(1, 1, 0))
#' pid <- ctrl_pid(kp = 2, ki = 0.5, kd = 1, derivative_filter = 10,
#'                 u_min = -1, u_max = 1, dt = 0.05)
#' sim <- ctrl_simulate_pid(plant, pid, reference = 1, steps = 400)
#' ctrl_performance(sim$time, sim$y, sim$r, u = sim$u)
#' plot(sim)
#' @export
ctrl_simulate_pid <- function(plant, pid, reference = 1, steps, disturbance = 0,
                              noise = 0, manual = NULL, x0 = NULL,
                              output = 1L) {
  if (!inherits(pid, "ctrl_pid")) ctrl_abort("`pid` must be a `ctrl_pid` object.")
  dt <- pid$dt
  if (dt <= 0) ctrl_abort("Simulation needs a PID controller with positive `dt`.")
  s <- as_ss(plant)
  if (ncol(s$B) != 1L) ctrl_abort("The PID controller drives a single input.")
  sd <- if (s$dt == 0) ctrl_c2d(s, dt) else s
  if (abs(sd$dt - dt) > 1e-12 * dt) {
    ctrl_abort("The discrete plant and controller sample intervals differ.")
  }
  Dp <- sd$D[output, 1]
  n <- nrow(sd$A)
  times <- (0:(steps - 1L)) * dt
  sig <- function(v, name) {
    if (is.function(v)) return(vapply(times, v, 0))
    v <- as.numeric(v)
    if (length(v) == 1L) return(rep(v, steps))
    if (length(v) != steps) ctrl_abort("`", name, "` must have length 1 or ", steps, ".")
    v
  }
  r <- sig(reference, "reference")
  d <- sig(disturbance, "disturbance")
  nz <- sig(noise, "noise")
  man <- if (is.null(manual)) rep(NA_real_, steps) else sig(manual, "manual")
  x <- if (is.null(x0)) rep(0, n) else as.numeric(x0)
  out <- data.frame(time = times, r = r, y = 0, y_true = 0, u = 0, u_unsat = 0,
                    P = 0, I = 0, D = 0, saturated = FALSE, mode = "auto",
                    stringsAsFactors = FALSE)
  st <- NULL
  for (k in seq_len(steps)) {
    mk <- if (is.na(man[k])) NULL else man[k]
    # y = C x + Dp (u + d): with feedthrough, u and y are solved together
    y0 <- drop(sd$C[output, , drop = FALSE] %*% x) + Dp * d[k]
    u_k <- if (Dp == 0) 0 else if (!is.null(mk)) {
      min(max(mk, pid$u_min), pid$u_max)
    } else pid_loop_solve(pid, st, r[k], y0 + nz[k], Dp)
    ytrue <- y0 + Dp * u_k
    ym <- ytrue + nz[k]
    step <- ctrl_pid_update(pid, st, r[k], ym, manual = mk)
    st <- step$state
    u <- step$u
    out$y[k] <- ym
    out$y_true[k] <- ytrue
    out$u[k] <- u
    out$u_unsat[k] <- step$u_unsat
    out$P[k] <- step$P
    out$I[k] <- step$I
    out$D[k] <- step$D
    out$saturated[k] <- step$saturated
    out$mode[k] <- if (is.na(man[k])) "auto" else "manual"
    x <- drop(sd$A %*% x + sd$B %*% (u + d[k]))
  }
  class(out) <- c("ctrl_pid_simulation", "data.frame")
  out
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_pid_simulation <- function(x, ...) {
  df <- as.data.frame(unclass(x))
  long <- rbind(
    data.frame(time = df$time, value = df$y, signal = "output y", panel = "output"),
    data.frame(time = df$time, value = df$r, signal = "reference r", panel = "output"),
    data.frame(time = df$time, value = df$u, signal = "control u", panel = "control"),
    data.frame(time = df$time, value = df$u_unsat, signal = "unsaturated",
               panel = "control")
  )
  ggplot2::ggplot(long, ggplot2::aes(x = .data$time, y = .data$value,
                                     color = .data$signal)) +
    ggplot2::geom_step() +
    ggplot2::facet_grid(panel ~ ., scales = "free_y") +
    ggplot2::labs(x = "time", y = NULL, color = NULL,
                  title = "PID closed-loop simulation") +
    ggplot2::theme_bw()
}

#' PID tuning rules
#'
#' * `"ziegler_nichols"`: the classical frequency-response rule from the
#'   ultimate gain \eqn{K_u} (the gain margin of the plant) and the ultimate
#'   period \eqn{P_u = 2\pi/\omega_{180}}: \eqn{k_p = 0.6 K_u},
#'   \eqn{T_i = P_u/2}, \eqn{T_d = P_u/8} (PID) or \eqn{k_p = 0.45 K_u},
#'   \eqn{T_i = P_u / 1.2} (PI). It gives aggressive, poorly damped loops
#'   and is included as a reference point.
#' * `"simc"`: Skogestad's SIMC rules for a first-order plus dead-time model
#'   \eqn{K e^{-\theta s}/(\tau s + 1)} with closed-loop time constant
#'   \eqn{\tau_c} (default \eqn{\theta}): \eqn{k_p = \tau / (K(\tau_c +
#'   \theta))}, \eqn{T_i = \min(\tau, 4(\tau_c + \theta))} (PI).
#' * `"amigo"`: the AMIGO PI rule of Astrom and Hagglund for the same model.
#'
#' Supply `plant` for Ziegler-Nichols, or `fopdt = c(K, tau, theta)` (for
#' example from [ctrl_fit_fopdt()]) for SIMC and AMIGO. The rules are
#' heuristics; verify the tuned loop by simulation and margins.
#'
#' @param method Tuning rule.
#' @param plant Plant model (Ziegler-Nichols).
#' @param fopdt Numeric vector `c(K, tau, theta)`.
#' @param type `"PID"` or `"PI"` (Ziegler-Nichols only; SIMC and AMIGO
#'   give PI controllers).
#' @param tau_c SIMC closed-loop time constant.
#' @param ... Further arguments for [ctrl_pid()] (for example `dt`,
#'   `u_min`, `u_max`, `derivative_filter`).
#' @return A `ctrl_pid` object.
#' @examples
#' plant <- ctrl_tf(1, c(1, 3, 3, 1))
#' ctrl_pid_tune("ziegler_nichols", plant = plant, derivative_filter = 10)
#' ctrl_pid_tune("simc", fopdt = c(1, 2, 0.5))
#' @export
ctrl_pid_tune <- function(method = c("ziegler_nichols", "simc", "amigo"),
                          plant = NULL, fopdt = NULL, type = c("PID", "PI"),
                          tau_c = NULL, ...) {
  method <- match.arg(method)
  type <- match.arg(type)
  if (method == "ziegler_nichols") {
    if (is.null(plant)) ctrl_abort("Ziegler-Nichols tuning needs `plant`.")
    mg <- ctrl_margin(plant)
    if (!is.finite(mg$gain_margin) || is.na(mg$phase_crossover) ||
        mg$phase_crossover <= 0) {
      ctrl_abort("The plant has no finite ultimate gain (no phase crossover).")
    }
    Ku <- mg$gain_margin
    Pu <- 2 * pi / mg$phase_crossover
    if (type == "PID") {
      return(ctrl_pid(kp = 0.6 * Ku, ki = Pu / 2, kd = Pu / 8,
                      form = "standard", ...))
    }
    return(ctrl_pid(kp = 0.45 * Ku, ki = Pu / 1.2, form = "standard", ...))
  }
  if (is.null(fopdt) || length(fopdt) != 3L) {
    ctrl_abort("`fopdt` must be c(K, tau, theta).")
  }
  K <- fopdt[1]
  tau <- fopdt[2]
  theta <- fopdt[3]
  if (tau <= 0 || theta < 0 || K == 0) ctrl_abort("Invalid FOPDT parameters.")
  if (method == "simc") {
    tau_c <- tau_c %||% max(theta, 1e-3 * tau)
    kp <- tau / (K * (tau_c + theta))
    Ti <- min(tau, 4 * (tau_c + theta))
    return(ctrl_pid(kp = kp, ki = Ti, form = "standard", ...))
  }
  kp <- (0.15 + (0.35 - theta * tau / (theta + tau)^2) * tau / max(theta, 1e-12)) / K
  Ti <- 0.35 * theta + 13 * theta * tau^2 / (tau^2 + 12 * theta * tau + 7 * theta^2)
  ctrl_pid(kp = kp, ki = max(Ti, 1e-12), form = "standard", ...)
}

#' Fit a first-order plus dead-time model to a step response
#'
#' Estimate \eqn{K e^{-\theta s}/(\tau s + 1)} from a unit step response
#' with the two-point method (28.3 percent and 63.2 percent of the final
#' change). The response must be monotone enough for the crossing times to
#' be meaningful.
#'
#' @param response A `ctrl_response` from [ctrl_step()], or a model (its
#'   step response is computed).
#' @return A named vector `c(K, tau, theta)`.
#' @examples
#' ctrl_fit_fopdt(ctrl_tf(1, c(1, 3, 3, 1)))
#' @export
ctrl_fit_fopdt <- function(response) {
  if (is_ctrl_system(response)) response <- ctrl_step(response)
  t <- response$time
  y <- response$y[, 1, 1]
  K <- y[length(y)] - y[1]
  if (K == 0) ctrl_abort("The step response does not change.")
  frac <- (y - y[1]) / K
  cross <- function(level) {
    k <- which(frac >= level)[1]
    if (is.na(k) || k == 1L) return(t[max(1L, k)])
    t[k - 1L] + (level - frac[k - 1L]) / (frac[k] - frac[k - 1L]) * (t[k] - t[k - 1L])
  }
  t1 <- cross(0.283)
  t2 <- cross(0.632)
  tau <- 1.5 * (t2 - t1)
  theta <- max(0, t2 - tau)
  c(K = K, tau = tau, theta = theta)
}
