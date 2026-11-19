#' System norms
#'
#' * `"2"`: the H2 norm. For an asymptotically stable continuous model it is
#'   \eqn{\sqrt{\mathrm{tr}(C W_c C^\top)}} when `D = 0` and infinite
#'   otherwise; for a discrete model it is
#'   \eqn{\sqrt{\mathrm{tr}(C W_c C^\top + D D^\top)}}.
#' * `"inf"`: the H-infinity norm \eqn{\sup_\omega \sigma_{\max}(G(j\omega))}
#'   of an asymptotically stable model, computed with the
#'   Bruinsma-Steinbuch two-step algorithm. The lower bound comes from
#'   evaluations at candidate frequencies; each iteration tests a level
#'   \eqn{\gamma} by checking the Hamiltonian matrix for imaginary-axis
#'   eigenvalues, which locate the frequency intervals where
#'   \eqn{\sigma_{\max} > \gamma}. The iteration stops when the relative gap
#'   between the bounds is below `tol`. Discrete models are mapped to
#'   continuous time with the norm-preserving bilinear transformation
#'   \eqn{s = (z - 1)/(z + 1)}.
#'
#' Unstable models have infinite H2 and H-infinity norms; `Inf` is returned.
#'
#' @param sys A model.
#' @param type `"2"` or `"inf"`.
#' @param tol Relative accuracy for the H-infinity norm.
#' @return A number; for `"inf"` with attribute `frequency` (the peak
#'   frequency in rad per time unit) and `bounds`.
#' @references Bruinsma, N. A. and Steinbuch, M. (1990). A fast algorithm to
#'   compute the H-infinity norm of a transfer function matrix. *Systems and
#'   Control Letters*, 14(4), 287-293.
#' @examples
#' G <- ctrl_tf(1, c(1, 0.2, 1))
#' ctrl_norm(G, "inf")   # 1 / (2 * 0.1 * sqrt(1 - 0.01))
#' ctrl_norm(ctrl_tf(1, c(1, 1)), "2")  # sqrt(1/2)
#' @export
ctrl_norm <- function(sys, type = c("2", "inf"), tol = 1e-10) {
  type <- match.arg(type)
  s <- as_ss(sys)
  if (!ctrl_is_stable(s)) return(Inf)
  if (type == "2") {
    if (!nrow(s$A)) {
      return(if (s$dt > 0) sqrt(sum(s$D^2)) else if (all(s$D == 0)) 0 else Inf)
    }
    Wc <- ctrl_gram(s, "c")
    if (s$dt > 0) {
      return(sqrt(max(0, sum(diag(s$C %*% Wc %*% t(s$C) + s$D %*% t(s$D))))))
    }
    if (any(s$D != 0)) return(Inf)
    return(sqrt(max(0, sum(diag(s$C %*% Wc %*% t(s$C))))))
  }
  hinf_norm(prescale_ss(s), tol)
}

hinf_norm <- function(s, tol) {
  dt <- s$dt
  if (!nrow(s$A)) {
    out <- max(svd(s$D, 0, 0)$d, 0)
    attr(out, "frequency") <- 0
    attr(out, "bounds") <- c(out, out)
    return(out)
  }
  if (dt > 0) {
    # bilinear map with T = 2: s = (z - 1) / (z + 1), frequencies tan(theta/2)
    sc <- s
    sc$dt <- 2
    cs <- ctrl_d2c(sc, method = "tustin")
    res <- hinf_norm(cs, tol)
    wc <- attr(res, "frequency")
    attr(res, "frequency") <- if (is.finite(wc)) 2 * atan(wc) / dt else pi / dt
    return(res)
  }
  A <- s$A
  B <- s$B
  C <- s$C
  D <- s$D
  n <- nrow(A)
  smax <- function(w) {
    G <- eval_system(s, 1i * w)[, , 1]
    max(svd(as.matrix(G), 0, 0)$d)
  }
  p <- as.complex(eigen(A, only.values = TRUE)$values)
  cand <- c(0, Mod(p), abs(Im(p)))
  # a log grid over the pole range guards against candidates that all hit
  # transmission zeros (a lower bound of zero would end the iteration)
  mp <- Mod(p)[Mod(p) > 0]
  lo <- if (length(mp)) log10(min(mp)) - 2 else -2
  hi <- if (length(mp)) log10(max(mp)) + 2 else 2
  cand <- c(cand, 10^seq(lo, hi, length.out = 41) * 1.0173)
  cand <- unique(cand[is.finite(cand)])
  vals <- vapply(cand, smax, 0)
  dinf <- max(svd(D, 0, 0)$d, 0)
  glb <- max(vals, dinf)
  wpk <- if (max(vals) >= dinf) cand[which.max(vals)] else Inf
  if (glb == 0) {
    # zero only if the transfer matrix vanishes identically: all Markov
    # parameters C A^k B (k < n) are zero (A scaled to avoid overflow)
    a <- max(norm(A, "2"), 1)
    M <- B
    for (k in seq_len(n)) {
      if (max(abs(C %*% M)) > 1e-14 * max(norm(C, "2"), 1) * max(norm(B, "2"), 1)) {
        ctrl_abort("Internal error: the H-infinity norm search found no ",
                   "nonzero response although the model is not zero.")
      }
      M <- (A / a) %*% M
    }
    out <- 0
    attr(out, "frequency") <- 0
    attr(out, "bounds") <- c(0, 0)
    return(out)
  }
  gub <- Inf
  for (iter in seq_len(60L)) {
    gam <- (1 + 2 * tol) * glb
    R <- gam^2 * diag(ncol(B)) - t(D) %*% D
    Sm <- gam^2 * diag(nrow(C)) - D %*% t(D)
    Ri <- solve(R)
    Ah <- A + B %*% Ri %*% t(D) %*% C
    H <- rbind(cbind(Ah, gam * B %*% Ri %*% t(B)),
               cbind(-gam * t(C) %*% solve(Sm) %*% C, -t(Ah)))
    ev <- eigen(H, only.values = TRUE)$values
    scale <- max(1, max(Mod(ev)))
    imag_ev <- ev[abs(Re(ev)) <= 1e-8 * scale & Im(ev) >= 0]
    if (!length(imag_ev)) {
      gub <- gam
      break
    }
    ws <- sort(Im(imag_ev))
    if (length(ws) == 1L) {
      mids <- ws
    } else {
      mids <- (ws[-1] + ws[-length(ws)]) / 2
      mids <- c(mids, ws)
    }
    v <- vapply(mids, smax, 0)
    if (max(v) <= glb) {
      gub <- gam
      break
    }
    glb <- max(v)
    wpk <- mids[which.max(v)]
  }
  out <- (glb + min(gub, (1 + 2 * tol) * glb)) / 2
  attr(out, "frequency") <- wpk
  attr(out, "bounds") <- c(glb, gub)
  out
}

#' Sensitivity functions
#'
#' For a loop transfer function `L` closed with unity negative feedback,
#' \eqn{S = (I + L)^{-1}} and \eqn{T = L (I + L)^{-1} = I - S}. With a plant
#' `P` and controller `C` (loop \eqn{L = P C}), `ctrl_gangof4()` also
#' returns \eqn{S P = (I + PC)^{-1} P} (input disturbance to output) and
#' \eqn{C S} (output noise to control signal).
#'
#' @param L Loop transfer function.
#' @param P,C Plant and controller.
#' @return A list of models `S` and `T` (and `PS`, `CS`).
#' @examples
#' L <- ctrl_tf(10, c(1, 3, 3, 1))
#' sens <- ctrl_sensitivity(L)
#' ctrl_norm(sens$S, "inf")  # peak sensitivity Ms
#' @export
ctrl_sensitivity <- function(L) {
  L <- as_system(L)
  p <- if (is_ctrl_tf(L)) nrow(L$num) else nrow(as_ss(L)$C)
  I <- if (is_ctrl_tf(L) && p == 1L) ctrl_tf(1, 1, dt = L$dt) else
    ctrl_gain(diag(p), dt = L$dt)
  S <- ctrl_feedback(I, L)
  T <- ctrl_feedback(L, I)
  list(S = S, T = T)
}

#' @rdname ctrl_sensitivity
#' @export
ctrl_gangof4 <- function(P, C) {
  P <- as_system(P)
  C <- as_system(C)
  L <- ctrl_series(C, P)
  sens <- ctrl_sensitivity(L)
  # S is the output sensitivity (I + PC)^-1; the load disturbance map is
  # S P (signal flows through P first) and the noise-to-control map is C S
  list(S = sens$S, T = sens$T, PS = ctrl_series(P, sens$S),
       CS = ctrl_series(sens$S, C))
}

#' Disk margins
#'
#' The disk margin of a SISO loop `L` (unity negative feedback) is the
#' largest \eqn{\alpha} such that the loop stays stable for every complex
#' gain perturbation \eqn{f} in the disk
#' \deqn{D(\alpha, \sigma) = \left\{\frac{1 + \frac{1-\sigma}{2}\delta}
#' {1 - \frac{1+\sigma}{2}\delta} : |\delta| \le \alpha\right\},}
#' where the skew \eqn{\sigma} shifts the disk (\eqn{\sigma = 0}: balanced,
#' \eqn{\sigma = 1}: gain increases only up to the classical multiplicative
#' form). It equals \eqn{\alpha_{\max} = 1 / \|S + \frac{\sigma - 1}{2}\|_\infty}
#' and guarantees *simultaneous* gain and phase variations inside the disk,
#' a stronger statement than the classical margins. The equivalent gain
#' margin interval and phase margin are reported. The nominal closed loop
#' must be stable.
#'
#' @param L SISO loop transfer function.
#' @param skew Skew \eqn{\sigma}.
#' @return A list of class `ctrl_disk_margin` with `alpha`, `gain_margin`
#'   (lower and upper factors), `gain_margin_db`, `phase_margin` (degrees)
#'   and `frequency` (the critical frequency).
#' @references Seiler, P., Packard, A. and Gahinet, P. (2020). An
#'   introduction to disk margins. *IEEE Control Systems Magazine*, 40(5),
#'   78-95.
#' @examples
#' L <- ctrl_tf(25, c(1, 10, 10, 10))
#' ctrl_disk_margin(L)
#' @export
ctrl_disk_margin <- function(L, skew = 0) {
  L_orig <- as_system(L)
  L <- siso_loop_tf(L)
  check_scalar_number(skew, "skew")
  S <- ctrl_sensitivity(L)$S
  hidden_ok <- tryCatch(ctrl_is_stable(ctrl_feedback(L_orig, 1)),
                        error = function(e) TRUE)
  if (!ctrl_is_stable(S) || !hidden_ok) {
    ctrl_abort("The nominal closed loop is unstable; disk margins are zero.")
  }
  Ss <- ctrl_tf2ss(S)
  Ss$D <- Ss$D + (skew - 1) / 2
  nrm <- ctrl_norm(Ss, "inf")
  alpha <- 1 / as.numeric(nrm)
  gm <- disk_gain_phase(alpha, skew)
  lo <- gm$lo
  hi <- gm$hi
  pm <- gm$pm
  structure(list(alpha = alpha, skew = skew, gain_margin = c(lo, hi),
                 gain_margin_db = 20 * log10(c(lo, hi)), phase_margin = pm,
                 frequency = attr(nrm, "frequency")),
            class = "ctrl_disk_margin")
}

# Gain and phase ranges of the disk {(1 + a d) / (1 - b d) : |d| <= alpha},
# a = (1 - skew) / 2, b = (1 + skew) / 2. The image is a disk (alpha < 1/b)
# or the exterior of a disk (alpha > 1/b) whose real-axis diameter runs
# between f(-alpha) and f(alpha).
disk_gain_phase <- function(alpha, skew) {
  # image of the disk |delta| <= alpha under f(delta) = (1 + a delta) /
  # (1 - b delta); f(alpha) and f(-alpha) are the real boundary points
  a <- (1 - skew) / 2
  b <- (1 + skew) / 2
  if (!is.finite(alpha)) return(list(lo = 0, hi = Inf, pm = 180))
  pa <- abs(b) * alpha
  if (abs(pa - 1) < 1e-12) {
    # pole of f on the boundary: the image is a half plane through the
    # finite real boundary point
    if (b > 0) {
      e1 <- (1 - a * alpha) / (1 + b * alpha)
      return(list(lo = max(e1, 0), hi = Inf, pm = acos(max(-1, min(1, e1))) * 180 / pi))
    }
    return(list(lo = 0, hi = (1 + a * alpha) / (1 - b * alpha), pm = 180))
  }
  e <- c((1 - a * alpha) / (1 + b * alpha), (1 + a * alpha) / (1 - b * alpha))
  m1 <- min(e)
  m2 <- max(e)
  # pole of f outside the disk: the image is the disk with diameter
  # [m1, m2]; pole inside: the exterior of that disk
  bounded <- pa < 1
  if (bounded) {
    lo <- max(m1, 0)
    hi <- m2
  } else if (m2 < 1) {
    lo <- max(m2, 0)
    hi <- Inf
  } else {
    lo <- 0
    hi <- m1
  }
  cc <- (m1 + m2) / 2
  r <- (m2 - m1) / 2
  k <- 1 + cc^2 - r^2
  # unit circle points inside the region satisfy 2 c cos(phi) >= k (bounded)
  # or 2 c cos(phi) <= k (exterior)
  if (abs(cc) < 1e-14) {
    pm <- 180
  } else if ((bounded && cc > 0) || (!bounded && cc < 0)) {
    pm <- acos(max(-1, min(1, k / (2 * cc)))) * 180 / pi
  } else {
    pm <- 180
  }
  list(lo = lo, hi = hi, pm = pm)
}

#' @export
print.ctrl_disk_margin <- function(x, digits = 4, ...) {
  cat("Disk margin (skew ", x$skew, "): alpha = ", signif(x$alpha, digits),
      " at w = ", signif(x$frequency, digits), "\n", sep = "")
  cat("  gain margin: [", signif(x$gain_margin[1], digits), ", ",
      signif(x$gain_margin[2], digits), "] (", signif(x$gain_margin_db[1],
                                                       digits), " to ",
      signif(x$gain_margin_db[2], digits), " dB)\n", sep = "")
  cat("  phase margin: +/- ", signif(x$phase_margin, digits), " deg\n",
      sep = "")
  invisible(x)
}
