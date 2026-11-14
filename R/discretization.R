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
#'   linearly interpolated between samples.
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
           D = s$D + s$C %*% f$G2)
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
  near0 <- function(v) any(Mod(v) < 1e-12)
  s0 <- if (near0(zpk$poles) || near0(zpk$zeros)) 1i / dt else 0
  gc <- zpk$gain * prod(s0 - zpk$zeros) / prod(s0 - zpk$poles)
  z0 <- exp(s0 * dt)
  gd0 <- prod(z0 - zz) / prod(z0 - pz)
  k <- Re(gc / gd0)
  out <- ctrl_zpk2tf(ctrl_zpk(zz, pz, k, dt = dt))
  out$input_names <- sys$input_names
  out$output_names <- sys$output_names
  out
}

#' Convert a discrete time model to continuous time
#'
#' `"zoh"` inverts the zero-order-hold map with a matrix logarithm of
#' \eqn{\begin{pmatrix} A_d & B_d \\ 0 & I\end{pmatrix}}; it fails when
#' \eqn{A_d} has eigenvalues on the closed negative real axis (no real
#' logarithm exists, the discrete model is not the sampled version of any
#' real continuous model). `"tustin"` inverts the bilinear map and fails
#' when \eqn{A_d} has an eigenvalue at -1. Both are exact inverses of the
#' corresponding [ctrl_c2d()] methods when they apply.
#'
#' @param sys A discrete time model.
#' @param method `"zoh"` or `"tustin"`.
#' @return A continuous time model.
#' @examples
#' sysd <- ctrl_c2d(ctrl_ss(-1, 1, 1), 0.1)
#' ctrl_d2c(sysd)$A  # -1
#' @export
ctrl_d2c <- function(sys, method = c("zoh", "tustin")) {
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
    if (any(Mod(Im(ev)) < 1e-12 & Re(ev) <= 1e-12)) {
      ctrl_abort("Zero-order-hold inversion is impossible: A has an ",
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
