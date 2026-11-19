#' Discretize a continuous time model
#'
#' Convert a continuous time model to discrete time with sample interval
#' `dt`.
#'
#' * `"zoh"` (default): exact zero-order hold. The input is held constant
#'   between samples and
#'   \eqn{A_d = e^{A T}}, \eqn{B_d = \int_0^T e^{A\tau} B\,d\tau}, both
#'   read from the matrix exponential of the block matrix
#'   \eqn{\begin{pmatrix} A & B \\ 0 & 0 \end{pmatrix} T}. This is valid for
#'   singular `A`; the formula \eqn{A^{-1}(A_d - I)B} is never used.
#' * `"foh"`: exact first-order (triangle) hold, for inputs that are
#'   linearly interpolated between samples. The discrete state is
#'   \eqn{\xi_k = x_k - \Gamma_2 u_k} (not the physical state), so a
#'   physical initial state \eqn{x_0} corresponds to
#'   \eqn{\xi_0 = x_0 - \Gamma_2 u_0}; the matrix \eqn{\Gamma_2} is stored
#'   in the `foh_state_map` attribute of the result.
#' * `"tustin"`: bilinear transformation \eqn{s = \frac{2}{T}\frac{z-1}{z+1}},
#'   optionally prewarped so that the frequency response matches exactly at
#'   `prewarp` rad per time unit.
#' * `"matched"`: matched pole-zero mapping \eqn{z = e^{sT}} for SISO models,
#'   with zeros at infinity mapped to \eqn{z = -1} (all but one, so the
#'   result stays strictly proper when the original is) and the gain
#'   matched at DC (or at `s = j/T` when there is a pole or zero at 0).
#'
#' Matrix exponentials use [expm::expm()] (scaling and squaring with Pade
#' approximation).
#'
#' @param sys A continuous time model.
#' @param dt Positive sample interval.
#' @param method One of `"zoh"`, `"foh"`, `"tustin"` or `"matched"`.
#' @param prewarp Optional prewarping frequency for `"tustin"`.
#' @return A discrete time model of the same class as `sys` (`"matched"`
#'   returns a transfer function).
#' @examples
#' sys <- ctrl_ss(-1, 1, 1)
#' sysd <- ctrl_c2d(sys, 0.1)
#' sysd$A  # exp(-0.1)
#' ctrl_c2d(ctrl_tf(1, c(1, 1, 0)), 0.5, method = "tustin")
#' @export
ctrl_c2d <- function(sys, dt, method = c("zoh", "foh", "tustin", "matched"),
                     prewarp = NULL) {
  method <- match.arg(method)
  assert_system(sys)
  check_scalar_number(dt, "dt", positive = TRUE)
  if (sys$dt > 0) ctrl_abort("`sys` is already a discrete time model.")
  if (method == "matched") return(c2d_matched(as_tf(sys), dt))
  was_tf <- is_ctrl_tf(sys) || inherits(sys, "ctrl_zpk")
  s <- as_ss(sys)
  n <- nrow(s$A)
  m <- ncol(s$B)
  if (n == 0L) {
    out <- ctrl_gain(s$D, dt = dt, input_names = s$input_names,
                     output_names = s$output_names)
    return(if (was_tf) ctrl_ss2tf(out) else out)
  }
  mats <- switch(method,
    zoh = zoh_matrices(s$A, s$B, dt),
    foh = {
      f <- foh_matrices(s$A, s$B, dt)
      list(A = f$Phi, B = f$G1 - f$G2 + f$Phi %*% f$G2, C = s$C,
           D = s$D + s$C %*% f$G2, G2 = f$G2)
    },
    tustin = {
      if (!is.null(prewarp)) {
        check_scalar_number(prewarp, "prewarp", positive = TRUE)
        if (prewarp * dt / 2 >= pi / 2) {
          ctrl_abort("`prewarp` must be below the Nyquist frequency pi / dt.")
        }
      }
      h <- if (is.null(prewarp)) dt else 2 / prewarp * tan(prewarp * dt / 2)
      ima <- diag(n) - 0.5 * h * s$A
      rc <- rcond(ima)
      if (rc < 1e-12) {
        ctrl_abort("Tustin discretization is undefined: A has an eigenvalue at ",
                   "2 / dt (rcond ", signif(rc, 3), ").")
      }
      Ad <- solve(ima, diag(n) + 0.5 * h * s$A)
      Bd <- solve(ima, h * s$B)
      Cd <- t(solve(t(ima), t(s$C)))
      list(A = Ad, B = Bd, C = Cd, D = s$D + 0.5 * s$C %*% Bd)
    }
  )
  out <- ctrl_ss(mats$A, mats$B, mats$C %||% s$C, mats$D %||% s$D, dt = dt,
                 state_names = s$state_names, input_names = s$input_names,
                 output_names = s$output_names)
  attr(out, "discretization") <- list(method = method, prewarp = prewarp)
  if (method == "foh") attr(out, "foh_state_map") <- mats$G2
  if (was_tf) ctrl_ss2tf(out) else out
}

zoh_matrices <- function(A, B, h) {
  n <- nrow(A)
  m <- ncol(B)
  M <- matrix(0, n + m, n + m)
  M[seq_len(n), seq_len(n)] <- A * h
  M[seq_len(n), n + seq_len(m)] <- B * h
  E <- expm::expm(M)
  list(A = E[seq_len(n), seq_len(n), drop = FALSE],
       B = E[seq_len(n), n + seq_len(m), drop = FALSE])
}

# Blocks of expm([A h, B h, 0; 0, 0, I; 0, 0, 0]). With linear interpolation
# of the input between samples, x(h) = Phi x0 + G1 u0 + G2 (u1 - u0).
foh_matrices <- function(A, B, h) {
  n <- nrow(A)
  m <- ncol(B)
  M <- matrix(0, n + 2L * m, n + 2L * m)
  M[seq_len(n), seq_len(n)] <- A * h
  M[seq_len(n), n + seq_len(m)] <- B * h
  M[n + seq_len(m), n + m + seq_len(m)] <- diag(m)
  E <- expm::expm(M)
  list(Phi = E[seq_len(n), seq_len(n), drop = FALSE],
       G1 = E[seq_len(n), n + seq_len(m), drop = FALSE],
       G2 = E[seq_len(n), n + m + seq_len(m), drop = FALSE])
}

c2d_matched <- function(sys, dt) {
  if (length(sys$num) != 1L) {
    ctrl_abort("Matched pole-zero discretization is only available for SISO ",
               "models.")
  }
  zpk <- ctrl_tf2zpk(sys)
  pz <- exp(zpk$poles * dt)
  zz <- exp(zpk$zeros * dt)
  reldeg <- length(zpk$poles) - length(zpk$zeros)
  if (reldeg > 1L) zz <- c(zz, rep(-1 + 0i, reldeg - 1L))
  # match the low-frequency gain after removing s^k and (z - 1)^k factors:
  # near s = 0, z - 1 ~ s dt
  # origin roots are classified once, relative to the root scale, and the
  # same classification is used on both sides; 1 - exp(p dt) is computed
  # with expm1 to avoid cancellation for slow roots
  rs <- Mod(c(zpk$zeros, zpk$poles))
  rsc <- if (length(rs) && max(rs) > 0) max(rs) else 1
  at0 <- function(v) Mod(v) <= 1e-12 * rsc
  oz <- at0(zpk$zeros)
  op <- at0(zpk$poles)
  kz <- sum(oz)
  kp <- sum(op)
  gc <- zpk$gain * prod(-zpk$zeros[!oz]) / prod(-zpk$poles[!op])
  # exp(z) - 1 for complex z without cancellation
  cexpm1 <- function(z) {
    x <- Re(z)
    y <- Im(z)
    complex(real = expm1(x) * cos(y) - 2 * sin(y / 2)^2, imaginary = exp(x) * sin(y))
  }
  gd <- prod(-cexpm1(zpk$zeros[!oz] * dt)) * 2^max(0L, reldeg - 1L) /
    prod(-cexpm1(zpk$poles[!op] * dt))
  k <- Re(gc / (gd * dt^(kz - kp)))
  out <- ctrl_zpk2tf(ctrl_zpk(zz, pz, k, dt = dt))
  out$input_names <- sys$input_names
  out$output_names <- sys$output_names
  out
}

#' Convert a discrete time model to continuous time
#'
#' `"zoh"` inverts the zero-order-hold map with a matrix logarithm of
#' \eqn{\begin{pmatrix} A_d & B_d \\ 0 & I\end{pmatrix}}; it fails when
#' \eqn{A_d} has eigenvalues on the closed negative real axis. Such
#' matrices may still have real logarithms (for example from pairs of
#' Jordan blocks), but the principal logarithm is not real and other
#' branches are not supported. `"tustin"` inverts the bilinear map and fails
#' when \eqn{A_d} has an eigenvalue at -1. Both are exact inverses of the
#' corresponding [ctrl_c2d()] methods when they apply.
#'
#' @param sys A discrete time model.
#' @param method `"zoh"` or `"tustin"`.
#' @param prewarp Prewarping frequency used in [ctrl_c2d()] (`"tustin"` only).
#' @return A continuous time model.
#' @examples
#' sysd <- ctrl_c2d(ctrl_ss(-1, 1, 1), 0.1)
#' ctrl_d2c(sysd)$A  # -1
#' @export
ctrl_d2c <- function(sys, method = c("zoh", "tustin"), prewarp = NULL) {
  method <- match.arg(method)
  assert_system(sys)
  if (sys$dt == 0) ctrl_abort("`sys` is already a continuous time model.")
  was_tf <- is_ctrl_tf(sys) || inherits(sys, "ctrl_zpk")
  s <- as_ss(sys)
  T <- s$dt
  n <- nrow(s$A)
  m <- ncol(s$B)
  if (n == 0L) {
    out <- ctrl_gain(s$D, input_names = s$input_names,
                     output_names = s$output_names)
    return(if (was_tf) ctrl_ss2tf(out) else out)
  }
  if (method == "zoh") {
    ev <- eigen(s$A, only.values = TRUE)$values
    # the principal logarithm needs no eigenvalue on the closed negative
    # real axis; small positive eigenvalues (fast stable poles) are fine
    if (any(abs(Im(ev)) <= 1e-12 * pmax(Mod(ev), 1e-300) & Re(ev) <= 0)) {
      ctrl_abort("Zero-order-hold inversion is not supported: A has an ",
                 "eigenvalue on the closed negative real axis.")
    }
    M <- rbind(cbind(s$A, s$B), cbind(matrix(0, m, n), diag(m)))
    L <- expm::logm(M) / T
    if (is.complex(L)) L <- Re(L)
    A <- L[seq_len(n), seq_len(n), drop = FALSE]
    B <- L[seq_len(n), n + seq_len(m), drop = FALSE]
    C <- s$C
    D <- s$D
  } else {
    ipa <- s$A + diag(n)
    if (rcond(ipa) < 1e-12) {
      ctrl_abort("Tustin inversion is impossible: A has an eigenvalue at -1.")
    }
    if (!is.null(prewarp)) {
      check_scalar_number(prewarp, "prewarp", positive = TRUE)
      if (prewarp * T / 2 >= pi / 2) {
        ctrl_abort("`prewarp` must be below the Nyquist frequency pi / dt.")
      }
      T <- 2 / prewarp * tan(prewarp * T / 2)
    }
    iinv <- solve(ipa)
    A <- (2 / T) * (s$A - diag(n)) %*% iinv
    B <- 2 * iinv %*% s$B / T
    C <- 2 * s$C %*% iinv
    D <- s$D - 0.5 * C %*% s$B
  }
  out <- ctrl_ss(A, B, C, D, state_names = s$state_names,
                 input_names = s$input_names, output_names = s$output_names)
  if (was_tf) ctrl_ss2tf(out) else out
}
