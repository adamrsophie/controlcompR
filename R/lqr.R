#' Linear quadratic regulator
#'
#' Design the state feedback \eqn{u = -K x} that minimizes, over inputs
#' that stabilize the closed loop,
#' \deqn{J = \int_0^\infty \left(x^\top Q x + u^\top R u + 2 x^\top N u
#' \right) dt}
#' for a continuous time model, or
#' \deqn{J = \sum_{k=0}^\infty \left(x_k^\top Q x_k + u_k^\top R u_k +
#' 2 x_k^\top N u_k\right)}
#' for a discrete time model. The gain is
#' \eqn{K = R^{-1}(B^\top P + N^\top)} (continuous) or
#' \eqn{K = (R + B^\top P B)^{-1}(B^\top P A + N^\top)} (discrete), where
#' `P` is the stabilizing solution of [ctrl_care()] or [ctrl_dare()].
#'
#' **Sign convention.** The control law is \eqn{u = -K x}, so the closed
#' loop is \eqn{A - B K}.
#'
#' Prerequisites are checked: `R` must be positive definite, the joint
#' weight \eqn{[Q, N; N^\top, R]} positive semidefinite and \eqn{(A, B)}
#' stabilizable. Without stabilization as a requirement the minimum can be
#' smaller: an unstable mode that the cost does not see (for example with
#' `Q = 0`) is cheapest left alone. The two coincide when the cost is
#' detectable. The optimality statement applies to the linear,
#' unconstrained, exactly known model; LQR is not optimal for nonlinear,
#' constrained or uncertain plants, and the guaranteed margins of
#' continuous time LQR (infinite gain margin, 60 degree phase margin per
#' loop when \eqn{N = 0} and `R` is diagonal) hold only with full state
#' measurement.
#'
#' @param sys A model.
#' @param Q State weight (`n x n`). A scalar or vector is expanded to a
#'   diagonal matrix.
#' @param R Input weight (`m x m`). A scalar or vector is expanded to a
#'   diagonal matrix.
#' @param N Optional cross weight (`n x m`).
#' @return An object of class `c("ctrl_lqr", "ctrl_controller")` with `K`,
#'   `P`, `poles` (closed-loop eigenvalues), `residual` (Riccati residual),
#'   the weights and the plant.
#' @seealso [ctrl_closed_loop()], [ctrl_dlqr_finite()], [ctrl_lqg()]
#' @examples
#' plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
#' controller <- ctrl_lqr(plant, Q = diag(c(10, 1)), R = 1)
#' controller$K
#' closed <- ctrl_closed_loop(plant, controller)
#' ctrl_stability(closed)
#' @export
ctrl_lqr <- function(sys, Q, R, N = NULL) {
  s <- as_ss(sys)
  n <- nrow(s$A)
  m <- ncol(s$B)
  if (!n) ctrl_abort("LQR needs a model with states.")
  Q <- expand_weight(Q, n, "Q")
  R <- expand_weight(R, m, "R")
  if (!is.null(N)) N <- as_real_matrix(N, "N", n, m)
  ctr <- ctrl_controllability(s)
  if (!ctr$stabilizable) {
    ctrl_abort("(A, B) is not stabilizable; uncontrollable unstable modes: ",
               paste(format(signif(ctr$uncontrollable_modes, 4)),
                     collapse = ", "), ".")
  }
  ric <- if (s$dt > 0) ctrl_dare(s$A, s$B, Q, R, N) else
    ctrl_care(s$A, s$B, Q, R, N)
  K <- ric$K
  dimnames(K) <- list(s$input_names, s$state_names)
  structure(
    list(K = K, P = ric$P, poles = ric$poles, residual = ric$residual,
         Q = Q, R = R, N = N %||% matrix(0, n, m), plant = s, dt = s$dt,
         riccati = ric),
    class = c("ctrl_lqr", "ctrl_controller")
  )
}

expand_weight <- function(W, k, name) {
  if (is.null(W)) ctrl_abort("`", name, "` is required.")
  if (is.null(dim(W))) {
    if (length(W) == 1L) W <- diag(W, k)
    else if (length(W) == k) W <- diag(W, k)
    else ctrl_abort("`", name, "` must be a scalar, a length ", k,
                    " vector or a ", k, " x ", k, " matrix.")
  }
  W <- as_real_matrix(W, name)
  check_dims(W, name, k, k)
  check_symmetric(W, name)
}

#' @export
print.ctrl_lqr <- function(x, digits = 4, ...) {
  cat("LQR controller (", time_domain_label(x$dt), "), u = -K x\n", sep = "")
  cat("K =\n")
  print(signif(x$K, digits))
  cat("closed-loop poles: ", paste(format(signif(x$poles, digits)),
                                   collapse = ", "), "\n", sep = "")
  cat("Riccati residual: ", signif(x$residual, 3), "\n", sep = "")
  invisible(x)
}

#' Finite-horizon discrete LQR
#'
#' Solve the finite-horizon problem
#' \deqn{\min \sum_{k=0}^{N-1} (x_k^\top Q x_k + u_k^\top R u_k) +
#' x_N^\top P_N x_N}
#' with the backward Riccati recursion
#' \deqn{P_k = Q + A^\top P_{k+1} A - A^\top P_{k+1} B (R + B^\top P_{k+1}
#' B)^{-1} B^\top P_{k+1} A.}
#' The optimal input is the time-varying feedback \eqn{u_k = -K_k x_k}. This
#' is the exact reference for unconstrained linear MPC with the same
#' weights and horizon.
#'
#' @param sys A discrete time model.
#' @param Q,R Weights (scalars and vectors are expanded to diagonals).
#' @param horizon Number of steps `N`.
#' @param P_terminal Terminal weight \eqn{P_N} (defaults to `Q`).
#' @return A list with `K` (list of `N` gains, `K[[1]]` applies at `k = 0`)
#'   and `P` (list of `N + 1` cost matrices, `P[[1]]` at `k = 0`).
#' @examples
#' sys <- ctrl_ss(1, 1, 1, dt = 1)
#' fin <- ctrl_dlqr_finite(sys, 1, 1, horizon = 30)
#' fin$K[[1]]  # approaches the infinite-horizon gain
#' ctrl_lqr(sys, 1, 1)$K
#' @export
ctrl_dlqr_finite <- function(sys, Q, R, horizon, P_terminal = NULL) {
  s <- as_ss(sys)
  if (s$dt == 0) ctrl_abort("`sys` must be a discrete time model.")
  n <- nrow(s$A)
  m <- ncol(s$B)
  Q <- check_psd(expand_weight(Q, n, "Q"), "Q")
  R <- expand_weight(R, m, "R")
  check_pd(R, "R")
  if (!is.numeric(horizon) || length(horizon) != 1L || horizon < 1 ||
      horizon != round(horizon)) {
    ctrl_abort("`horizon` must be a positive integer.")
  }
  PN <- if (is.null(P_terminal)) Q else
    check_psd(expand_weight(P_terminal, n, "P_terminal"), "P_terminal")
  P <- vector("list", horizon + 1L)
  K <- vector("list", horizon)
  P[[horizon + 1L]] <- PN
  A <- s$A
  B <- s$B
  for (k in horizon:1) {
    Pn <- P[[k + 1L]]
    G <- R + t(B) %*% Pn %*% B
    K[[k]] <- solve(G, t(B) %*% Pn %*% A)
    P[[k]] <- symmetrize(Q + t(A) %*% Pn %*% (A - B %*% K[[k]]))
  }
  list(K = K, P = P)
}

#' Closed-loop systems
#'
#' Build the closed loop of a plant with a controller.
#'
#' * With an LQR or other static state feedback (`ctrl_lqr`, or a matrix
#'   `K`), the control law is \eqn{u = -K x + v}, where `v` is an external
#'   input. The closed loop is
#'   \eqn{\dot x = (A - BK) x + B v}, \eqn{y = (C - DK) x + D v}, with
#'   outputs `y` followed (when `include_control = TRUE`) by the control
#'   signal `u`.
#' * With an LQG controller (see [ctrl_lqg()]) the controller is a dynamic
#'   system from the measurement to the input, connected in positive
#'   feedback `u = K_c(y)` (the controller already contains the minus
#'   sign), with an additive input disturbance `v` at the plant input.
#'   States are the plant states followed by the controller states.
#' * With a dynamic `ctrl_ss` or `ctrl_tf` controller `Cs` acting on the
#'   error \eqn{r - y}, the loop is the unity negative feedback
#'   \eqn{y = (I + G C_s)^{-1} G C_s r}.
#'
#' @param plant Plant model.
#' @param controller A `ctrl_lqr` object, a gain matrix `K`, a `ctrl_lqg`
#'   object, or a `ctrl_ss` / `ctrl_tf` error-feedback controller.
#' @param include_control Append the control signal to the outputs.
#' @return A `ctrl_ss` model.
#' @examples
#' plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
#' cl <- ctrl_closed_loop(plant, ctrl_lqr(plant, diag(2), 1))
#' ctrl_poles(cl)
#' @export
ctrl_closed_loop <- function(plant, controller, include_control = TRUE) {
  s <- as_ss(plant)
  if (inherits(controller, "ctrl_lqg")) {
    return(lqg_closed_loop(s, controller, include_control))
  }
  if (is_ctrl_system(controller)) {
    return(ctrl_feedback(ctrl_series(controller, s), ctrl_gain(diag(nrow(s$C)),
                                                               dt = s$dt)))
  }
  K <- if (inherits(controller, "ctrl_controller")) controller$K else controller
  K <- as_real_matrix(K, "K", ncol(s$B), nrow(s$A))
  check_dims(K, "K", ncol(s$B), nrow(s$A))
  if (inherits(controller, "ctrl_controller") && !is.null(controller$dt) &&
      controller$dt != s$dt) {
    ctrl_abort("The controller and plant have different time domains.")
  }
  Acl <- s$A - s$B %*% K
  Ccl <- s$C - s$D %*% K
  Dcl <- s$D
  out_names <- s$output_names
  if (include_control) {
    Ccl <- rbind(Ccl, -K)
    Dcl <- rbind(Dcl, diag(ncol(s$B)))
    out_names <- c(out_names, paste0("ctrl_", s$input_names))
  }
  ctrl_ss(Acl, s$B, Ccl, Dcl, dt = s$dt, state_names = s$state_names,
          input_names = paste0("v_", s$input_names), output_names = out_names)
}

#' Reference (feedforward) gain for state feedback
#'
#' For the law \eqn{u = -K x + N_r r}, compute \eqn{N_r} so that the closed
#' loop has unit DC gain from `r` to the outputs `y` (square case, least
#' squares otherwise). This gives zero steady-state error for constant
#' references only for the exact nominal model; use integral action or
#' offset-free MPC for robustness to model error.
#'
#' @param plant Plant model.
#' @param K Gain matrix or `ctrl_lqr` object.
#' @return The `m x p` matrix \eqn{N_r}.
#' @examples
#' plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
#' lqr <- ctrl_lqr(plant, diag(2), 1)
#' ctrl_reference_gain(plant, lqr)
#' @export
ctrl_reference_gain <- function(plant, K) {
  cl <- ctrl_closed_loop(plant, K, include_control = FALSE)
  G0 <- ctrl_dcgain(cl)
  if (nrow(G0) == ncol(G0)) {
    safe_solve(G0, diag(nrow(G0)), "closed-loop DC gain")
  } else {
    MASS_ginv(G0)
  }
}

MASS_ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X)
  pos <- s$d > max(tol * s$d[1L], 0)
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}
