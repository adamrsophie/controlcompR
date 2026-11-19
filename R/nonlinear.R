#' Nonlinear state-space models
#'
#' Define a nonlinear model
#' \eqn{\dot x = f(t, x, u, \theta)} (continuous time, `dt = NULL`) or
#' \eqn{x_{k+1} = f(k, x_k, u_k, \theta)} (discrete time, `dt > 0`) with
#' output \eqn{y = h(t, x, u, \theta)}.
#'
#' @param update Function `function(t, x, u, params)` returning the state
#'   derivative (continuous) or the next state (discrete).
#' @param output Function `function(t, x, u, params)` returning the output;
#'   `NULL` outputs the full state.
#' @param nstates,ninputs Dimensions.
#' @param noutputs Number of outputs (defaults to `nstates` without
#'   `output`; otherwise evaluated once).
#' @param dt `NULL` for continuous time or a positive sample interval.
#' @param params List of parameters passed to `update` and `output`.
#' @param state_names,input_names,output_names Optional names.
#' @return An object of class `ctrl_nlsys`.
#' @seealso [ctrl_nl_simulate()], [ctrl_linearize()], [ctrl_equilibrium()]
#' @examples
#' pendulum <- ctrl_nlsys(
#'   function(t, x, u, p) c(x[2], -p$g / p$l * sin(x[1]) - p$c * x[2] + u),
#'   output = function(t, x, u, p) x[1],
#'   nstates = 2, ninputs = 1, params = list(g = 9.81, l = 1, c = 0.5)
#' )
#' pendulum
#' @export
ctrl_nlsys <- function(update, output = NULL, nstates, ninputs, noutputs = NULL,
                       dt = NULL, params = list(), state_names = NULL,
                       input_names = NULL, output_names = NULL) {
  if (!is.function(update)) ctrl_abort("`update` must be a function.")
  if (!is.null(output) && !is.function(output)) ctrl_abort("`output` must be a function.")
  n <- as_count(nstates, "nstates")
  m <- if (identical(ninputs, 0) || identical(ninputs, 0L)) 0L else as_count(ninputs, "ninputs")
  dt <- if (is.null(dt)) 0 else match_time_domain(dt)
  if (is.null(output)) {
    output <- function(t, x, u, params) x
    p <- n
  } else {
    p <- noutputs %||% length(output(0, rep(0, n), rep(0, m), params))
  }
  fx <- update(0, rep(0, n), rep(0, m), params)
  if (!is.numeric(fx) || length(fx) != n) {
    ctrl_abort("`update` must return a numeric vector of length ", n, ".")
  }
  structure(list(update = update, output = output, n = n, m = m, p = p, dt = dt,
                 params = params,
                 state_names = state_names %||% default_names("x", n),
                 input_names = input_names %||% default_names("u", m),
                 output_names = output_names %||% default_names("y", p)),
            class = "ctrl_nlsys")
}

#' @export
print.ctrl_nlsys <- function(x, ...) {
  cat("Nonlinear model (", time_domain_label(x$dt), "): ", x$n, " states, ",
      x$m, " inputs, ", x$p, " outputs\n", sep = "")
  invisible(x)
}

nl_f <- function(sys, t, x, u) as.numeric(sys$update(t, x, u, sys$params))
nl_h <- function(sys, t, x, u) as.numeric(sys$output(t, x, u, sys$params))

#' Simulate a nonlinear model
#'
#' Continuous time models are integrated with `deSolve` (default method
#' `"lsoda"`, an adaptive stiff/nonstiff solver); discrete models are
#' iterated. The input can be
#'
#' * a vector or matrix with one row per time point, held constant between
#'   points (zero-order hold, integrated interval by interval);
#' * a function `input(t)` evaluated continuously by the integrator;
#' * a feedback law `controller(t, x)` evaluated at the time points and held
#'   in between, which simulates a sampled-data controller.
#'
#' @param sys A `ctrl_nlsys` model.
#' @param x0 Initial state.
#' @param times Time points (for discrete models, sample indices are
#'   `times / dt`).
#' @param input Input (see above) or `NULL` for zero input.
#' @param controller Optional feedback law `function(t, x)`.
#' @param method `deSolve` integration method.
#' @param rtol,atol Integration tolerances.
#' @return An object of class `ctrl_nl_response` with `time`, `x`, `u` and
#'   `y` matrices.
#' @examples
#' pendulum <- ctrl_nlsys(
#'   function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
#'   nstates = 2, ninputs = 1
#' )
#' sim <- ctrl_nl_simulate(pendulum, x0 = c(2, 0), times = seq(0, 10, 0.05))
#' plot(sim)
#' @export
ctrl_nl_simulate <- function(sys, x0, times, input = NULL, controller = NULL,
                             method = "lsoda", rtol = 1e-8, atol = 1e-10) {
  if (!inherits(sys, "ctrl_nlsys")) ctrl_abort("`sys` must be a `ctrl_nlsys` model.")
  x0 <- as.numeric(x0)
  if (length(x0) != sys$n) ctrl_abort("`x0` must have length ", sys$n, ".")
  times <- as.numeric(times)
  Tn <- length(times)
  if (Tn < 2L || any(diff(times) <= 0)) ctrl_abort("`times` must be increasing.")
  m <- sys$m
  X <- matrix(0, Tn, sys$n)
  Um <- matrix(0, Tn, m)
  X[1, ] <- x0
  uf <- NULL
  if (is.function(input)) {
    uf <- input
  } else if (!is.null(input)) {
    Um <- input_matrix(input, times, m)
  }
  if (sys$dt > 0) {
    x <- x0
    for (k in seq_len(Tn)) {
      X[k, ] <- x
      if (!is.null(controller)) Um[k, ] <- as.numeric(controller(times[k], x))
      if (!is.null(uf)) Um[k, ] <- as.numeric(uf(times[k]))
      if (k < Tn) {
        ratio <- (times[k + 1L] - times[k]) / sys$dt
        steps <- round(ratio)
        if (abs(ratio - steps) > 1e-8 * max(1, ratio)) {
          ctrl_abort("`times` must be spaced by integer multiples of the ",
                     "sample interval dt = ", sys$dt, ".")
        }
        for (j in seq_len(steps)) x <- nl_f(sys, times[k] / sys$dt + j - 1, x, Um[k, ])
      }
    }
  } else if (!is.null(uf) && is.null(controller)) {
    out <- deSolve::ode(x0, times, function(t, x, p) list(nl_f(sys, t, x, as.numeric(uf(t)))),
                        NULL, method = method, rtol = rtol, atol = atol)
    X <- unname(out[, -1, drop = FALSE])
    for (k in seq_len(Tn)) Um[k, ] <- as.numeric(uf(times[k]))
  } else {
    x <- x0
    for (k in seq_len(Tn)) {
      X[k, ] <- x
      if (!is.null(controller)) Um[k, ] <- as.numeric(controller(times[k], x))
      if (k < Tn) {
        uk <- Um[k, ]
        out <- deSolve::ode(x, times[k:(k + 1L)],
                            function(t, x, p) list(nl_f(sys, t, x, uk)), NULL,
                            method = method, rtol = rtol, atol = atol)
        if (nrow(out) < 2L || anyNA(out[2, ])) {
          ctrl_abort("Integration failed on [", times[k], ", ", times[k + 1L], "].")
        }
        x <- unname(out[2, -1])
      }
    }
  }
  # discrete models receive the sample index, continuous models the time
  tk <- if (sys$dt > 0) times / sys$dt else times
  Y <- t(vapply(seq_len(Tn), function(k) nl_h(sys, tk[k], X[k, ], Um[k, ]),
                numeric(sys$p)))
  if (sys$p == 1L) Y <- matrix(Y, ncol = 1)
  colnames(X) <- sys$state_names
  colnames(Y) <- sys$output_names
  if (m) colnames(Um) <- sys$input_names
  structure(list(time = times, x = X, u = Um, y = Y, dt = sys$dt),
            class = "ctrl_nl_response")
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_nl_response <- function(x, ...) {
  rows <- list()
  for (j in seq_len(ncol(x$x))) {
    rows[[length(rows) + 1L]] <- data.frame(time = x$time, value = x$x[, j],
                                            signal = colnames(x$x)[j], panel = "states")
  }
  for (j in seq_len(ncol(x$u))) {
    rows[[length(rows) + 1L]] <- data.frame(time = x$time, value = x$u[, j],
                                            signal = colnames(x$u)[j], panel = "inputs")
  }
  ggplot2::ggplot(do.call(rbind, rows),
                  ggplot2::aes(x = .data$time, y = .data$value, color = .data$signal)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(panel ~ ., scales = "free_y") +
    ggplot2::labs(x = time_label(x$dt), y = NULL, color = NULL,
                  title = "Nonlinear simulation") +
    ggplot2::theme_bw()
}

num_jacobian <- function(fun, x, method = "central") {
  f0 <- fun(x)
  J <- matrix(0, length(f0), length(x))
  for (i in seq_along(x)) {
    if (method == "central") {
      h <- .Machine$double.eps^(1 / 3) * max(1, abs(x[i]))
      e <- replace(numeric(length(x)), i, h)
      J[, i] <- (fun(x + e) - fun(x - e)) / (2 * h)
    } else {
      h <- sqrt(.Machine$double.eps) * max(1, abs(x[i]))
      e <- replace(numeric(length(x)), i, h)
      J[, i] <- (fun(x + e) - f0) / h
    }
  }
  J
}

#' Linearize a nonlinear model
#'
#' Compute the Jacobians \eqn{A = \partial f/\partial x},
#' \eqn{B = \partial f/\partial u}, \eqn{C = \partial h/\partial x},
#' \eqn{D = \partial h/\partial u} at an operating point. The result is a
#' local model in deviation variables \eqn{\delta x = x - x_0},
#' \eqn{\delta u = u - u_0}; it is exact only to first order and only near
#' the operating point. Jacobians are computed by central finite
#' differences (step \eqn{\epsilon^{1/3}\max(1, |x_i|)}, error
#' \eqn{O(h^2)}) unless analytic Jacobians are supplied. The method used is
#' recorded in the result.
#'
#' At an operating point that is not an equilibrium, the linearization has
#' an additional constant term \eqn{f(x_0, u_0)} (returned as `drift`) that
#' the `ctrl_ss` model does not include.
#'
#' @param sys A `ctrl_nlsys` model.
#' @param x0,u0 Operating point.
#' @param t Time at which to linearize.
#' @param method `"central"` or `"forward"` differences.
#' @param jacobian Optional function `function(t, x, u, params)` returning
#'   `list(A = , B = )` (and optionally `C`, `D`) analytically.
#' @return A `ctrl_ss` model with attributes `operating_point`, `drift` and
#'   `method`.
#' @examples
#' pendulum <- ctrl_nlsys(
#'   function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
#'   output = function(t, x, u, p) x[1], nstates = 2, ninputs = 1
#' )
#' lin <- ctrl_linearize(pendulum, x0 = c(pi, 0), u0 = 0)
#' ctrl_poles(lin)  # upright position: unstable
#' @export
ctrl_linearize <- function(sys, x0, u0 = numeric(sys$m), t = 0,
                           method = c("central", "forward"), jacobian = NULL) {
  method <- match.arg(method)
  x0 <- as.numeric(x0)
  u0 <- as.numeric(u0)
  if (length(x0) != sys$n || length(u0) != sys$m) {
    ctrl_abort("The operating point must have ", sys$n, " states and ", sys$m, " inputs.")
  }
  jac <- if (!is.null(jacobian)) jacobian(t, x0, u0, sys$params) else list()
  A <- jac$A %||% num_jacobian(function(x) nl_f(sys, t, x, u0), x0, method)
  B <- jac$B %||% if (sys$m) num_jacobian(function(u) nl_f(sys, t, x0, u), u0, method) else
    matrix(0, sys$n, 0)
  C <- jac$C %||% num_jacobian(function(x) nl_h(sys, t, x, u0), x0, method)
  D <- jac$D %||% if (sys$m) num_jacobian(function(u) nl_h(sys, t, x0, u), u0, method) else
    matrix(0, sys$p, 0)
  out <- ctrl_ss(as.matrix(A), as.matrix(B), as.matrix(C), as.matrix(D),
                 dt = if (sys$dt > 0) sys$dt else NULL,
                 state_names = sys$state_names, input_names = sys$input_names,
                 output_names = sys$output_names)
  drift <- nl_f(sys, t, x0, u0) - if (sys$dt > 0) x0 else 0
  attr(out, "operating_point") <- list(x = x0, u = u0, t = t)
  attr(out, "drift") <- drift
  attr(out, "method") <- if (is.null(jacobian)) paste(method, "finite differences") else
    "analytic (user supplied)"
  out
}

#' Equilibrium points
#'
#' Solve \eqn{f(x, u) = 0} (continuous) or \eqn{f(x, u) = x} (discrete)
#' with a damped Newton method using finite-difference Jacobians. By
#' default the inputs are fixed at `u0` and all states are unknown; with
#' `free_inputs`, the listed input components are also unknowns, and with
#' `fixed_states`, the listed states are held at their values in `x0` (for
#' example to find the input that holds a given operating point). When
#' there are more unknowns than equations, minimum-norm Newton steps are
#' used and the solution depends on the initial guess.
#' Local stability is assessed from the linearization: an equilibrium with
#' a linearization on the stability boundary is reported as
#' `"inconclusive"`, since nonlinear terms decide stability there.
#'
#' @param sys A `ctrl_nlsys` model.
#' @param x0 Initial guess for the state.
#' @param u0 Input (or initial guess for the free inputs).
#' @param free_inputs Indices of inputs to solve for.
#' @param fixed_states Indices of states held at `x0`.
#' @param tol Residual tolerance.
#' @param maxit Maximum Newton iterations.
#' @return A list with `x`, `u`, `residual`, `converged`, `iterations`,
#'   `eigenvalues` (of the linearization) and `stability`.
#' @examples
#' pendulum <- ctrl_nlsys(
#'   function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
#'   nstates = 2, ninputs = 1
#' )
#' ctrl_equilibrium(pendulum, x0 = c(3, 0))$x
#' @export
ctrl_equilibrium <- function(sys, x0, u0 = numeric(sys$m), free_inputs = integer(),
                             fixed_states = integer(), tol = 1e-10, maxit = 100L) {
  x <- as.numeric(x0)
  u <- as.numeric(u0)
  free <- as.integer(free_inputs)
  xfree <- setdiff(seq_len(sys$n), as.integer(fixed_states))
  nxf <- length(xfree)
  res_fun <- function(z) {
    xx <- x
    xx[xfree] <- z[seq_len(nxf)]
    uu <- u
    uu[free] <- z[nxf + seq_along(free)]
    nl_f(sys, 0, xx, uu) - if (sys$dt > 0) xx else 0
  }
  z <- c(x[xfree], u[free])
  r <- res_fun(z)
  it <- 0L
  while (sqrt(sum(r^2)) > tol && it < maxit) {
    J <- num_jacobian(res_fun, z)
    step <- tryCatch(-drop(MASS_ginv(J) %*% r), error = function(e) rep(0, length(z)))
    lam <- 1
    nr <- sqrt(sum(r^2))
    repeat {
      zn <- z + lam * step
      rn <- res_fun(zn)
      if (all(is.finite(rn)) && sqrt(sum(rn^2)) < (1 - 1e-4 * lam) * nr) break
      lam <- lam / 2
      if (lam < 1e-10) break
    }
    if (lam < 1e-10) break
    z <- zn
    r <- rn
    it <- it + 1L
  }
  x[xfree] <- z[seq_len(nxf)]
  u[free] <- z[nxf + seq_along(free)]
  lin <- ctrl_linearize(sys, x, u)
  ev <- ctrl_poles(lin)
  converged <- sqrt(sum(r^2)) <= tol
  marg <- if (sys$dt > 0) max(Mod(ev)) - 1 else max(Re(ev))
  bnd <- 1e-8 * max(1, max(Mod(ev)))
  # stability is a property of an equilibrium; without one, only the
  # Jacobian eigenvalues at the last iterate are reported
  stab <- if (!converged) "not assessed (no equilibrium found)" else
    if (marg < -bnd) "locally asymptotically stable" else if (marg > bnd)
    "unstable" else "inconclusive (linearization on the stability boundary)"
  list(x = x, u = u, residual = sqrt(sum(r^2)), converged = converged,
       iterations = it, eigenvalues = ev, stability = stab)
}

#' Discretize a continuous nonlinear model
#'
#' Build a discrete time `ctrl_nlsys` whose update integrates the
#' continuous model over one sample interval with the input held constant
#' (`deSolve`, tight tolerances, or a fixed number of classical
#' Runge-Kutta steps with `method = "rk4"`).
#'
#' @param sys Continuous `ctrl_nlsys`.
#' @param dt Sample interval.
#' @param method `"ode"` (deSolve lsoda) or `"rk4"`.
#' @param substeps RK4 steps per sample.
#' @return A discrete `ctrl_nlsys`.
#' @examples
#' pendulum <- ctrl_nlsys(
#'   function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
#'   nstates = 2, ninputs = 1
#' )
#' pd <- ctrl_nl_discretize(pendulum, 0.05, method = "rk4")
#' @export
ctrl_nl_discretize <- function(sys, dt, method = c("ode", "rk4"), substeps = 4L) {
  method <- match.arg(method)
  if (sys$dt > 0) ctrl_abort("The model is already discrete.")
  check_scalar_number(dt, "dt", positive = TRUE)
  substeps <- as_count(substeps, "substeps")
  f <- sys$update
  out_fun <- sys$output
  step <- if (method == "rk4") {
    h <- dt / substeps
    function(k, x, u, params) {
      t <- k * dt
      for (s in seq_len(substeps)) {
        k1 <- f(t, x, u, params)
        k2 <- f(t + h / 2, x + h / 2 * k1, u, params)
        k3 <- f(t + h / 2, x + h / 2 * k2, u, params)
        k4 <- f(t + h, x + h * k3, u, params)
        x <- x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t <- t + h
      }
      x
    }
  } else {
    function(k, x, u, params) {
      o <- deSolve::ode(x, c(k * dt, (k + 1) * dt),
                        function(t, x, p) list(as.numeric(f(t, x, u, params))), NULL,
                        rtol = 1e-10, atol = 1e-12)
      unname(o[2, -1])
    }
  }
  ctrl_nlsys(step, output = function(k, x, u, params) out_fun(k * dt, x, u, params),
             nstates = sys$n, ninputs = sys$m, noutputs = sys$p, dt = dt,
             params = sys$params, state_names = sys$state_names,
             input_names = sys$input_names, output_names = sys$output_names)
}

nl_filter_inputs <- function(sys, y, u, Qn, Rn, x0, P0) {
  if (!inherits(sys, "ctrl_nlsys") || sys$dt == 0) {
    ctrl_abort("Nonlinear filters need a discrete time `ctrl_nlsys`; use ctrl_nl_discretize().")
  }
  Y <- measurement_matrix(y, sys$p)
  Tn <- nrow(Y)
  U <- if (is.null(u)) matrix(0, Tn, sys$m) else input_matrix(u, seq_len(Tn), sys$m)
  list(Y = Y, U = U, Tn = Tn,
       Qn = check_psd(expand_weight(Qn, sys$n, "Qn"), "Qn"),
       Rn = check_pd(expand_weight(Rn, sys$p, "Rn"), "Rn"),
       x = if (is.null(x0)) rep(0, sys$n) else as.numeric(x0),
       P = if (is.null(P0)) diag(sys$n) else
         check_psd(expand_weight(P0, sys$n, "P0"), "P0"))
}

#' Extended and unscented Kalman filters
#'
#' Nonlinear state estimation for a discrete time `ctrl_nlsys` with
#' additive noise, \eqn{x_{k+1} = f(k, x_k, u_k) + w_k},
#' \eqn{y_k = h(k, x_k, u_k) + v_k}, using the conventions of
#' [ctrl_kalman_filter()] (`x0`, `P0` are the prior at the first sample;
#' each step updates with \eqn{y_k} and then predicts; rows of `y` that are
#' entirely `NA` skip the update). Row `k` of the data is passed to `f` and
#' `h` as the sample index `k - 1`, the 0-based convention of
#' [ctrl_nl_simulate()].
#'
#' * `ctrl_ekf()` linearizes `f` and `h` at the current estimate with
#'   central finite differences and uses the Joseph-form covariance update.
#' * `ctrl_ukf()` propagates \eqn{2n + 1} sigma points (scaled unscented
#'   transform with parameters `alpha`, `beta`, `kappa`); it needs no
#'   Jacobians and is exact for linear models.
#'
#' Both are approximations for nonlinear models: their covariances are not
#' guaranteed to be consistent, and they can diverge for strong
#' nonlinearities or poor initialization. Check the normalized innovation
#' squared (`nis`), whose mean should be close to the number of outputs.
#'
#' @param sys Discrete time `ctrl_nlsys`.
#' @param y Measurements (`T x p`).
#' @param u Inputs (`T x m`) or `NULL`.
#' @param Qn,Rn Process and measurement noise covariances.
#' @param x0,P0 Prior mean and covariance.
#' @param alpha,beta,kappa Unscented transform parameters.
#' @return A list of class `c("ctrl_nl_filter", "ctrl_kalman_filter")` with
#'   `filtered`, `P_filtered`, `innovations`, `nis` and `loglik`. It has no
#'   predicted covariances, so [ctrl_kalman_smoother()] does not accept it.
#' @examples
#' f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] - 0.1 * sin(x[1]))
#' h <- function(k, x, u, p) x[1]
#' sys <- ctrl_nlsys(f, h, nstates = 2, ninputs = 0, dt = 0.1)
#' set.seed(1)
#' x <- c(1, 0)
#' y <- numeric(100)
#' for (k in 1:100) {
#'   y[k] <- x[1] + rnorm(1, sd = 0.05)
#'   x <- f(k, x, NULL, NULL) + rnorm(2, sd = 0.001)
#' }
#' ekf <- ctrl_ekf(sys, y, Qn = 1e-6, Rn = 0.0025, x0 = c(0.8, 0))
#' ukf <- ctrl_ukf(sys, y, Qn = 1e-6, Rn = 0.0025, x0 = c(0.8, 0))
#' tail(cbind(ekf$filtered, ukf$filtered), 2)
#' @export
ctrl_ekf <- function(sys, y, u = NULL, Qn, Rn, x0 = NULL, P0 = NULL) {
  z <- nl_filter_inputs(sys, y, u, Qn, Rn, x0, P0)
  n <- sys$n
  p <- sys$p
  x <- z$x
  P <- z$P
  # histories are written in place in this loop (no per-sample copies)
  Tn <- z$Tn
  xf <- matrix(0, Tn, n)
  Pf <- array(0, c(n, n, Tn))
  innov <- matrix(NA_real_, Tn, p)
  nis <- rep(NA_real_, Tn)
  loglik <- 0
  for (k in seq_len(Tn)) {
    tk <- k - 1L  # sample index passed to the model, 0-based as in simulation
    uk <- z$U[k, ]
    yk <- z$Y[k, ]
    if (all(is.finite(yk))) {
      H <- num_jacobian(function(xx) nl_h(sys, tk, xx, uk), x)
      e <- yk - nl_h(sys, tk, x, uk)
      S <- symmetrize(H %*% P %*% t(H) + z$Rn)
      K <- P %*% t(H) %*% solve(S)
      x <- x + drop(K %*% e)
      IKH <- diag(n) - K %*% H
      P <- symmetrize(IKH %*% P %*% t(IKH) + K %*% z$Rn %*% t(K))
      st <- nl_innovation_stats(e, S)
      innov[k, ] <- e
      nis[k] <- st[["nis"]]
      loglik <- loglik + st[["ll"]]
    }
    xf[k, ] <- x
    Pf[, , k] <- P
    Fj <- num_jacobian(function(xx) nl_f(sys, tk, xx, uk), x)
    x <- nl_f(sys, tk, x, uk)
    P <- symmetrize(Fj %*% P %*% t(Fj) + z$Qn)
  }
  nl_filter_finish(list(filtered = xf, P_filtered = Pf, innovations = innov,
                        nis = nis, loglik = loglik), sys)
}

#' @rdname ctrl_ekf
#' @export
ctrl_ukf <- function(sys, y, u = NULL, Qn, Rn, x0 = NULL, P0 = NULL, alpha = 1,
                     beta = 2, kappa = 0) {
  z <- nl_filter_inputs(sys, y, u, Qn, Rn, x0, P0)
  n <- sys$n
  p <- sys$p
  check_scalar_number(alpha, "alpha", positive = TRUE)
  check_scalar_number(beta, "beta")
  check_scalar_number(kappa, "kappa")
  lam <- alpha^2 * (n + kappa) - n
  if (n + lam <= 0) ctrl_abort("The unscented transform needs n + lambda > 0; increase `kappa`.")
  wm <- c(lam / (n + lam), rep(1 / (2 * (n + lam)), 2 * n))
  wc <- wm
  wc[1] <- wc[1] + (1 - alpha^2 + beta)
  sigma <- function(x, P) {
    Lc <- t(chol_psd(P))
    S <- sqrt(n + lam) * Lc
    cbind(x, x + S, x - S)
  }
  x <- z$x
  P <- z$P
  Tn <- z$Tn
  xf <- matrix(0, Tn, n)
  Pf <- array(0, c(n, n, Tn))
  innov <- matrix(NA_real_, Tn, p)
  nis <- rep(NA_real_, Tn)
  loglik <- 0
  for (k in seq_len(Tn)) {
    tk <- k - 1L
    uk <- z$U[k, ]
    if (all(is.finite(z$Y[k, ]))) {
      Xs <- sigma(x, P)
      Ys <- matrix(apply(Xs, 2, function(xx) nl_h(sys, tk, xx, uk)), nrow = p)
      yhat <- drop(Ys %*% wm)
      dY <- Ys - yhat
      dX <- Xs - x
      Pyy <- symmetrize(dY %*% (wc * t(dY)))
      S <- Pyy + z$Rn
      Pxy <- dX %*% (wc * t(dY))
      K <- Pxy %*% solve(S)
      e <- z$Y[k, ] - yhat
      x <- x + drop(K %*% e)
      # Joseph form with the statistically linearized output map
      # H = Pxy' P^+, plus the part of Pyy it does not explain (zero for
      # linear models). This equals P - K S K' algebraically but avoids its
      # cancellation when the posterior is much smaller than the prior.
      Hs <- t(Pxy) %*% MASS_ginv(P)
      IKH <- diag(n) - K %*% Hs
      Rx <- z$Rn + (Pyy - Hs %*% P %*% t(Hs))
      P <- symmetrize(IKH %*% P %*% t(IKH) + K %*% Rx %*% t(K))
      st <- nl_innovation_stats(e, S)
      innov[k, ] <- e
      nis[k] <- st[["nis"]]
      loglik <- loglik + st[["ll"]]
    }
    xf[k, ] <- x
    Pf[, , k] <- P
    Xs <- sigma(x, P)
    Fs <- matrix(apply(Xs, 2, function(xx) nl_f(sys, tk, xx, uk)), nrow = n)
    x <- drop(Fs %*% wm)
    dF <- Fs - x
    P <- symmetrize(dF %*% (wc * t(dF)) + z$Qn)
  }
  nl_filter_finish(list(filtered = xf, P_filtered = Pf, innovations = innov,
                        nis = nis, loglik = loglik), sys)
}

chol_psd <- function(P) {
  r <- tryCatch(chol(P), error = function(e) NULL)
  if (!is.null(r)) return(r)
  ev <- eigen(symmetrize(P), symmetric = TRUE)
  t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow(P)))
}

# Normalized innovation squared and log-likelihood term of one update.
nl_innovation_stats <- function(e, S) {
  cS <- chol(S)
  zz <- forwardsolve(t(cS), e)
  nis <- sum(zz^2)
  c(nis = nis, ll = -0.5 * (length(e) * log(2 * pi) + 2 * sum(log(diag(cS))) + nis))
}

nl_filter_finish <- function(out, sys) {
  colnames(out$filtered) <- sys$state_names
  out$plant <- sys
  class(out) <- c("ctrl_nl_filter", "ctrl_kalman_filter")
  out
}

#' Nonlinear model predictive control (experimental)
#'
#' Single-shooting nonlinear MPC for a discrete time `ctrl_nlsys`: at each
#' sample the input sequence minimizing
#' \deqn{\sum_{k=0}^{N-1} \|x_k - x^r\|_Q^2 + \|u_k - u^r\|_R^2 +
#' \|x_N - x^r\|_P^2 + \rho \sum_k \|\max(0, x_k - x_{\max},
#' x_{\min} - x_k)\|^2}
#' over the box \eqn{u_{\min} \le u_k \le u_{\max}} is computed with the
#' bounded quasi-Newton method L-BFGS-B (`stats::optim`). The gradient is
#' computed by a backward (adjoint) pass through the prediction, with the
#' stage Jacobians of the model approximated by central differences, so it
#' costs \eqn{O(N (n + m))} model evaluations. Input bounds are hard; state bounds are enforced only
#' through the quadratic penalty \eqn{\rho} and can be violated.
#'
#' **Limitations.** The problem is nonconvex: the solver returns a local
#' minimum that depends on the initial guess (the shifted previous
#' solution is used), and no stability or
#' recursive feasibility guarantee is provided. Each step reports the
#' optimizer convergence code; nonzero codes are treated as failures.
#'
#' @param sys Discrete time `ctrl_nlsys`.
#' @param horizon Prediction horizon.
#' @param Q,R,P Weights (`P` defaults to `Q`).
#' @param u_min,u_max Input bounds.
#' @param x_min,x_max State bounds (soft, penalized).
#' @param penalty Penalty weight for state bound violations.
#' @param maxit Maximum optimizer iterations per step.
#' @return An object of class `ctrl_nmpc`.
#' @examples
#' f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] + 0.1 * (-sin(x[1]) + u))
#' sys <- ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1)
#' nmpc <- ctrl_nmpc(sys, horizon = 15, Q = diag(c(10, 1)), R = 0.1,
#'                   u_min = -2, u_max = 2)
#' sim <- ctrl_simulate_nmpc(nmpc, x0 = c(1, 0), steps = 30)
#' tail(sim$x, 1)
#' @export
ctrl_nmpc <- function(sys, horizon, Q, R, P = NULL, u_min = -Inf, u_max = Inf,
                      x_min = -Inf, x_max = Inf, penalty = 1e4, maxit = 200L) {
  if (!inherits(sys, "ctrl_nlsys") || sys$dt == 0) {
    ctrl_abort("NMPC needs a discrete time `ctrl_nlsys`; use ctrl_nl_discretize().")
  }
  n <- sys$n
  m <- sys$m
  if (!m) ctrl_abort("NMPC needs a model with inputs.")
  Q <- check_psd(expand_weight(Q, n, "Q"), "Q")
  R <- check_pd(expand_weight(R, m, "R"), "R")
  P <- if (is.null(P)) Q else check_psd(expand_weight(P, n, "P"), "P")
  structure(list(plant = sys, horizon = as_count(horizon, "horizon"), Q = Q, R = R,
                 P = P, u_min = expand_bound(u_min, m, "u_min"),
                 u_max = expand_bound(u_max, m, "u_max"),
                 x_min = expand_bound(x_min, n, "x_min"),
                 x_max = expand_bound(x_max, n, "x_max"), penalty = penalty,
                 maxit = maxit, dt = sys$dt),
            class = c("ctrl_nmpc", "ctrl_controller"))
}

#' @rdname ctrl_nmpc
#' @param nmpc A `ctrl_nmpc` object.
#' @param x Current state.
#' @param x_ref,u_ref State and input references: vectors (constant) or
#'   matrices whose row `i` is the reference at sample `i - 1` (the last row
#'   is repeated beyond the end), or functions of the 1-based sample index.
#' @param u_init Initial guess (`N x m` matrix).
#' @param k Time index.
#' @export
ctrl_nmpc_step <- function(nmpc, x, x_ref = NULL, u_ref = NULL, u_init = NULL,
                           k = 0L) {
  s <- nmpc$plant
  n <- s$n
  m <- s$m
  N <- nmpc$horizon
  # references: vectors (constant) or matrices whose row i is sample i - 1
  # of the whole run (the window starts at row k + 1)
  Xr <- ref_rows(x_ref, k + 1L, n, N + 1L, "x_ref")
  Ur <- ref_rows(u_ref, k + 1L, m, N, "u_ref")
  u0 <- if (is.null(u_init)) {
    t(pmin(pmax(t(Ur), nmpc$u_min), nmpc$u_max))
  } else as.matrix(u_init)
  if (!identical(dim(u0), c(N, m))) ctrl_abort("`u_init` must be a ", N, " x ", m, " matrix.")
  rollout <- function(uv) {
    U <- matrix(uv, N, m)
    X <- matrix(0, N + 1L, n)
    X[1, ] <- x
    for (j in seq_len(N)) X[j + 1L, ] <- nl_f(s, k + j - 1L, X[j, ], U[j, ])
    list(X = X, U = U)
  }
  cost <- function(uv) {
    r <- rollout(uv)
    J <- 0
    for (j in seq_len(N)) {
      dx <- r$X[j, ] - Xr[j, ]
      du <- r$U[j, ] - Ur[j, ]
      J <- J + sum(dx * (nmpc$Q %*% dx)) + sum(du * (nmpc$R %*% du))
    }
    dN <- r$X[N + 1L, ] - Xr[N + 1L, ]
    J <- J + sum(dN * (nmpc$P %*% dN))
    v <- pmax(0, sweep(r$X[-1, , drop = FALSE], 2, nmpc$x_max),
              sweep(-r$X[-1, , drop = FALSE], 2, -nmpc$x_min))
    J + nmpc$penalty * sum(v^2)
  }
  # exact gradient of the single-shooting cost by a backward (adjoint)
  # pass; only the stage Jacobians are approximated, by central differences
  grad <- function(uv) {
    r <- rollout(uv)
    G <- matrix(0, N, m)
    pen_grad <- function(x) {
      2 * nmpc$penalty * (pmax(0, x - nmpc$x_max) - pmax(0, nmpc$x_min - x))
    }
    lam <- 2 * drop(nmpc$P %*% (r$X[N + 1L, ] - Xr[N + 1L, ])) + pen_grad(r$X[N + 1L, ])
    for (j in N:1) {
      xj <- r$X[j, ]
      uj <- r$U[j, ]
      tj <- k + j - 1L
      Fx <- num_jacobian(function(xx) nl_f(s, tj, xx, uj), xj)
      Fu <- num_jacobian(function(uu) nl_f(s, tj, xj, uu), uj)
      G[j, ] <- 2 * drop(nmpc$R %*% (uj - Ur[j, ])) + drop(t(Fu) %*% lam)
      lam <- 2 * drop(nmpc$Q %*% (xj - Xr[j, ])) + drop(t(Fx) %*% lam) +
        if (j > 1L) pen_grad(xj) else 0
    }
    as.vector(G)
  }
  lo <- rep(nmpc$u_min, each = N)
  hi <- rep(nmpc$u_max, each = N)
  start <- pmin(pmax(as.vector(u0), lo), hi)
  opt <- tryCatch(stats::optim(start, cost, grad, method = "L-BFGS-B", lower = lo,
                               upper = hi, control = list(maxit = nmpc$maxit)),
                  error = function(e) list(par = start, value = NA_real_,
                                           convergence = 99L, message = conditionMessage(e)))
  r <- rollout(opt$par)
  list(u = r$U[1, ], ok = opt$convergence == 0L && is.finite(opt$value),
       convergence = opt$convergence, message = opt$message, objective = opt$value,
       u_plan = r$U, x_plan = r$X)
}

#' @rdname ctrl_nmpc
#' @param x0 Initial state.
#' @param steps Number of closed-loop samples.
#' @param plant Simulation model (default: the prediction model). It must
#'   be a discrete time `ctrl_nlsys` with the controller's sample interval
#'   and dimensions; discretize a continuous plant with
#'   [ctrl_nl_discretize()] first.
#' @param on_failure `"stop"` or `"hold"` (apply the previous input). The
#'   held input is clipped to the input bounds; the log records the
#'   fallback, and `"hold (clipped)"` when clipping changed it. A held input
#'   does not guarantee that state constraints stay satisfied.
#' @param u0 Input held before the first sample, used by the `"hold"`
#'   fallback when the first optimization fails. The default is the
#'   admissible point of the input bounds closest to zero.
#' @export
ctrl_simulate_nmpc <- function(nmpc, x0, steps, x_ref = NULL, u_ref = NULL,
                               plant = NULL, on_failure = c("stop", "hold"),
                               u0 = NULL) {
  on_failure <- match.arg(on_failure)
  sim <- plant %||% nmpc$plant
  if (!inherits(sim, "ctrl_nlsys") || sim$dt == 0) {
    ctrl_abort("`plant` must be a discrete time `ctrl_nlsys`; use ctrl_nl_discretize().")
  }
  if (abs(sim$dt - nmpc$dt) > 1e-12 * nmpc$dt) {
    ctrl_abort("`plant` has sample interval ", sim$dt, " but the controller uses ",
               nmpc$dt, ".")
  }
  if (sim$n != nmpc$plant$n || sim$m != nmpc$plant$m) {
    ctrl_abort("`plant` must have ", nmpc$plant$n, " states and ", nmpc$plant$m,
               " inputs like the prediction model.")
  }
  steps <- as_count(steps, "steps")
  n <- sim$n
  m <- sim$m
  X <- matrix(0, steps + 1L, n)
  U <- matrix(0, steps, m)
  x0 <- as.numeric(x0)
  if (length(x0) != n || any(!is.finite(x0))) {
    ctrl_abort("`x0` must be a finite vector of length ", n, ".")
  }
  X[1, ] <- x0
  log <- data.frame(k = seq_len(steps), ok = FALSE, convergence = 0L,
                    objective = NA_real_, fallback = "none", stringsAsFactors = FALSE)
  guess <- NULL
  u_prev <- if (is.null(u0)) pmin(pmax(0, nmpc$u_min), nmpc$u_max) else {
    u0 <- as.numeric(u0)
    if (length(u0) != m || any(!is.finite(u0))) {
      ctrl_abort("`u0` must be a finite vector of length ", m, ".")
    }
    if (any(u0 < nmpc$u_min | u0 > nmpc$u_max)) {
      ctrl_abort("`u0` must satisfy the input bounds.")
    }
    u0
  }
  for (k in seq_len(steps)) {
    st <- ctrl_nmpc_step(nmpc, X[k, ], x_ref, u_ref, guess, k - 1L)
    log$ok[k] <- st$ok
    log$convergence[k] <- st$convergence
    log$objective[k] <- st$objective
    if (st$ok) {
      u <- st$u
      guess <- rbind(st$u_plan[-1, , drop = FALSE], st$u_plan[nmpc$horizon, , drop = FALSE])
    } else if (on_failure == "stop") {
      ctrl_abort("NMPC optimization failed at sample ", k, " (code ", st$convergence,
                 ": ", st$message %||% "", ").")
    } else {
      u <- pmin(pmax(u_prev, nmpc$u_min), nmpc$u_max)
      log$fallback[k] <- if (any(u != u_prev)) "hold (clipped)" else "hold"
    }
    U[k, ] <- u
    X[k + 1L, ] <- nl_f(sim, k - 1L, X[k, ], u)
    u_prev <- u
  }
  colnames(X) <- sim$state_names
  colnames(U) <- sim$input_names
  list(time = (0:steps) * sim$dt, x = X, u = U, log = log, failures = sum(!log$ok))
}

#' Nonlinear trajectory optimization
#'
#' Compute an open-loop input sequence \eqn{u_0, \ldots, u_{N-1}} for a
#' discrete time `ctrl_nlsys` that minimizes
#' \deqn{\sum_{k=0}^{N-1} \|x_k - x^r_k\|_Q^2 + \|u_k - u^r_k\|_R^2 +
#' \|x_N - x^r_N\|_P^2 + \rho \sum_k \|\max(0, x_k - x_{\max},
#' x_{\min} - x_k)\|^2}
#' from a given initial state, by single shooting with L-BFGS-B and the
#' adjoint gradient of [ctrl_nmpc()]. Input bounds are hard; state bounds are
#' penalized.
#'
#' **Terminal constraint.** A target `x_final` is imposed as the equality
#' \eqn{x_N = x_f} by an augmented Lagrangian method: each outer iteration
#' minimizes the cost plus
#' \eqn{\nu^\top (x_N - x_f) + \frac{\rho}{2}\|x_N - x_f\|^2} and then
#' updates \eqn{\nu \leftarrow \nu + \rho (x_N - x_f)} (multiplying
#' \eqn{\rho} by 10 when the violation shrinks by less than a factor 4),
#' until
#' \eqn{\|x_N - x_f\| \le} `terminal_tol`. A moderate \eqn{\rho}
#' (`terminal_weight`) keeps each subproblem well conditioned, unlike a
#' single large penalty. `P` is ignored when `x_final` is given.
#'
#' The problem is nonconvex and the result is a local minimum that depends
#' on `u_init`. Failures are reported in `converged` and `message`, never
#' hidden.
#'
#' @inheritParams ctrl_nmpc
#' @param x0 Initial state.
#' @param x_ref,u_ref References (see [ctrl_nmpc_step()]).
#' @param x_final Optional target state at the final sample.
#' @param terminal_weight Penalty parameter \eqn{\rho} for `x_final`.
#' @param terminal_tol Required accuracy of the terminal state.
#' @param max_outer Maximum number of multiplier updates.
#' @param u_init Initial guess (`N x m` matrix).
#' @return An object of class `ctrl_trajopt` with `time`, `x` (`N + 1` rows),
#'   `u` (`N` rows), `objective` (the stated cost, including the soft
#'   state bound penalty but not the terminal constraint terms), `cost`
#'   (its `quadratic` and `penalty` parts), `converged`, `convergence` (optimizer code of the last
#'   subproblem), `message`, `outer_iterations`, `multiplier` and
#'   `terminal_error` (distance to `x_final`, `NA` without one).
#' @examples
#' f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] + 0.1 * (-sin(x[1]) + u))
#' sys <- ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1)
#' tr <- ctrl_nl_trajopt(sys, x0 = c(0, 0), horizon = 40, Q = 0, R = 0.01,
#'                       x_final = c(pi / 2, 0), u_min = -2, u_max = 2)
#' tr$terminal_error
#' @export
ctrl_nl_trajopt <- function(sys, x0, horizon, Q, R, P = NULL, x_ref = NULL,
                            u_ref = NULL, x_final = NULL, terminal_weight = 100,
                            terminal_tol = 1e-6, max_outer = 30L,
                            u_min = -Inf, u_max = Inf, x_min = -Inf, x_max = Inf,
                            penalty = 1e4, maxit = 1000L, u_init = NULL) {
  N <- as_count(horizon, "horizon")
  n <- sys$n
  x0 <- as.numeric(x0)
  if (length(x0) != n || any(!is.finite(x0))) {
    ctrl_abort("`x0` must be a finite vector of length ", n, ".")
  }
  mk <- function(P) {
    ctrl_nmpc(sys, N, Q, R, P = P, u_min = u_min, u_max = u_max,
              x_min = x_min, x_max = x_max, penalty = penalty, maxit = maxit)
  }
  outer <- 0L
  nu <- NULL
  terr <- NA_real_
  if (is.null(x_final)) {
    st <- ctrl_nmpc_step(mk(P), x0, x_ref, u_ref, u_init, k = 0L)
    ok <- st$ok
  } else {
    x_final <- as.numeric(x_final)
    if (length(x_final) != n) ctrl_abort("`x_final` must have length ", n, ".")
    check_scalar_number(terminal_weight, "terminal_weight", positive = TRUE)
    rho <- terminal_weight
    nm <- mk(rho / 2 * diag(n))
    Xr <- ref_rows(x_ref, 1L, n, N + 1L, "x_ref")
    nu <- numeric(n)
    terr_prev <- Inf
    guess <- u_init
    ok <- FALSE
    for (outer in seq_len(as_count(max_outer, "max_outer"))) {
      # nu'(x - xf) + rho/2 |x - xf|^2 = rho/2 |x - xf + nu/rho|^2 + const
      Xr[N + 1L, ] <- x_final - nu / rho
      st <- ctrl_nmpc_step(nm, x0, Xr, u_ref, guess, k = 0L)
      guess <- st$u_plan
      dN <- st$x_plan[N + 1L, ] - x_final
      terr <- sqrt(sum(dN^2))
      if (terr <= terminal_tol) {
        ok <- st$ok
        break
      }
      nu <- nu + rho * dN
      # standard rule: raise rho when the violation did not shrink enough
      if (outer > 1L && terr > 0.25 * terr_prev) {
        rho <- 10 * rho
        nm <- mk(rho / 2 * diag(n))
      }
      terr_prev <- terr
    }
  }
  X <- st$x_plan
  U <- st$u_plan
  # report the cost of the original problem (without terminal terms)
  Xr0 <- ref_rows(x_ref, 1L, n, N + 1L, "x_ref")
  Ur0 <- ref_rows(u_ref, 1L, sys$m, N, "u_ref")
  Qm <- expand_weight(Q, n, "Q")
  Rm <- expand_weight(R, sys$m, "R")
  J <- 0
  for (j in seq_len(N)) {
    dx <- X[j, ] - Xr0[j, ]
    du <- U[j, ] - Ur0[j, ]
    J <- J + sum(dx * (Qm %*% dx)) + sum(du * (Rm %*% du))
  }
  if (is.null(x_final)) {
    Pm <- if (is.null(P)) Qm else expand_weight(P, n, "P")
    dN <- X[N + 1L, ] - Xr0[N + 1L, ]
    J <- J + sum(dN * (Pm %*% dN))
  }
  # the soft state bounds are part of the stated objective
  xmin <- expand_bound(x_min, n, "x_min")
  xmax <- expand_bound(x_max, n, "x_max")
  viol <- pmax(0, sweep(X[-1, , drop = FALSE], 2, xmax),
               sweep(-X[-1, , drop = FALSE], 2, -xmin))
  Jpen <- penalty * sum(viol^2)
  colnames(X) <- sys$state_names
  colnames(U) <- sys$input_names
  structure(list(time = (0:N) * sys$dt, x = X, u = U, objective = J + Jpen,
                 cost = c(quadratic = J, penalty = Jpen),
                 converged = ok, convergence = st$convergence,
                 message = if (!ok && !is.null(x_final) && !(terr <= terminal_tol))
                   "terminal constraint not met within max_outer iterations" else st$message,
                 outer_iterations = outer, multiplier = nu, terminal_error = terr),
            class = "ctrl_trajopt")
}

#' @export
print.ctrl_trajopt <- function(x, digits = 4, ...) {
  cat("Trajectory optimization: ", nrow(x$u), " steps, objective ",
      signif(x$objective, digits), ", converged: ", x$converged, "\n", sep = "")
  if (is.finite(x$terminal_error)) {
    cat("  terminal error: ", signif(x$terminal_error, digits), "\n", sep = "")
  }
  invisible(x)
}
