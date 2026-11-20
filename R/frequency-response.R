#' Frequency response
#'
#' Evaluate \eqn{G(j\omega)} for a continuous time model or
#' \eqn{G(e^{j\omega\,dt})} for a discrete time model. For discrete models
#' frequencies above the Nyquist frequency \eqn{\pi/dt} alias and produce a
#' warning.
#'
#' State-space models are balanced (diagonal scaling) and evaluated by
#' solving \eqn{(j\omega I - A) X = B} at each frequency (no matrix
#' inverse). Transfer functions are evaluated as ratios of their numerator
#' and denominator polynomials (Horner's rule), which can lose accuracy for
#' high order polynomials; convert to state space with [ctrl_tf2ss()] for
#' those.
#'
#' @param sys A model.
#' @param w Frequencies in rad per time unit. If `NULL`, a logarithmic grid
#'   is chosen from the pole and zero locations (see [ctrl_freq_grid()]).
#' @return A list of class `ctrl_freqresp` with `w` and `response`, a
#'   complex array of dimension outputs x inputs x frequencies.
#' @examples
#' G <- ctrl_tf(1, c(1, 1))
#' fr <- ctrl_freqresp(G, c(0.1, 1, 10))
#' Mod(fr$response[1, 1, ])
#' @export
ctrl_freqresp <- function(sys, w = NULL) {
  assert_system(sys)
  w <- w %||% ctrl_freq_grid(sys)
  if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0)) {
    ctrl_abort("`w` must be nonnegative finite frequencies.")
  }
  if (sys$dt > 0 && any(w > pi / sys$dt * (1 + 1e-12))) {
    ctrl_warn("Some frequencies exceed the Nyquist frequency pi / dt = ",
              signif(pi / sys$dt, 4), " and alias.")
  }
  z <- if (sys$dt > 0) exp(1i * w * sys$dt) else 1i * w
  resp <- eval_system(sys, z)
  structure(list(w = w, response = resp, dt = sys$dt,
                 input_names = sys$input_names,
                 output_names = sys$output_names),
            class = "ctrl_freqresp")
}

# Evaluate a model at complex points z (s or z plane); returns p x m x nz.
eval_system <- function(sys, z) {
  if (is_ctrl_tf(sys)) return(tf_eval(sys, z))
  if (inherits(sys, "ctrl_zpk")) return(tf_eval(ctrl_zpk2tf(sys), z))
  s <- prescale_ss(sys)
  n <- nrow(s$A)
  p <- nrow(s$C)
  m <- ncol(s$B)
  out <- array(0 + 0i, c(p, m, length(z)))
  if (!n) {
    for (k in seq_along(z)) out[, , k] <- s$D
    return(out)
  }
  # Hessenberg form would be faster for many frequencies; a dense solve per
  # frequency is O(n^3) each.
  I <- diag(n)
  for (k in seq_along(z)) {
    M <- z[k] * I - s$A
    # no rcond threshold: a badly scaled but nonsingular resolvent (a pole
    # at 1e100 next to poles at the origin) is solved accurately by LU;
    # only an exactly singular or overflowing solve gives NA
    X <- tryCatch(solve(M, s$B, tol = 0), error = function(e) NULL)
    Y <- if (is.null(X)) NULL else s$C %*% X + s$D
    out[, , k] <- if (is.null(Y) || any(!is.finite(Y))) NA_complex_ else Y
  }
  out
}

#' Frequency grid for a model
#'
#' A logarithmic grid from one to two decades below the smallest to one to
#' two decades above the largest nonzero pole and zero magnitude (the ends
#' are rounded outward to whole decades), refined near lightly damped poles. For
#' discrete models the grid ends at the Nyquist frequency.
#'
#' @param sys A model.
#' @param n Number of points in the base grid.
#' @return A sorted vector of positive frequencies.
#' @examples
#' range(ctrl_freq_grid(ctrl_tf(1, c(1, 0.2, 100))))
#' @export
ctrl_freq_grid <- function(sys, n = 400L) {
  p <- ctrl_poles(sys)
  z <- tryCatch(ctrl_zeros(sys), error = function(e) complex())
  pts <- c(p, z)
  if (sys$dt > 0) {
    pts <- pts[Mod(pts) > 1e-12]
    pts <- log(pts) / sys$dt
  }
  mags <- Mod(pts)
  mags <- mags[mags > 1e-12 & is.finite(mags)]
  lo <- if (length(mags)) floor(log10(min(mags))) - 1 else -2
  hi <- if (length(mags)) ceiling(log10(max(mags))) + 1 else 2
  if (sys$dt > 0) hi <- min(hi, log10(pi / sys$dt))
  if (hi <= lo) lo <- hi - 3
  w <- 10^seq(lo, hi, length.out = n)
  # refine around lightly damped resonances
  res <- pts[abs(Im(pts)) > 1e-12]
  for (r in res) {
    wn <- Mod(r)
    zeta <- -Re(r) / wn
    if (abs(zeta) < 0.5) {
      bw <- max(abs(zeta), 1e-4) * wn
      w <- c(w, wn + bw * seq(-5, 5, length.out = 41))
    }
  }
  w <- w[w > 0]
  if (sys$dt > 0) w <- w[w <= pi / sys$dt]
  sort(unique(w))
}

unwrap_phase <- function(ph) {
  if (length(ph) < 2L) return(ph)
  d <- diff(ph)
  adj <- -2 * pi * round(d / (2 * pi))
  ph + c(0, cumsum(adj))
}

#' Bode, Nyquist and Nichols data
#'
#' `ctrl_bode()` returns magnitude (absolute and dB) and phase (degrees)
#' over a frequency grid. Phase is unwrapped along the grid: the result is
#' continuous only when the grid resolves every rapid phase change, so
#' check lightly damped modes with a denser `w`. `ctrl_nyquist()` returns
#' the complex response for positive and negative frequencies and, for SISO
#' models without poles on the stability boundary, the winding number of
#' \eqn{1 + L} around the origin (the number of clockwise encirclements of
#' \eqn{-1}); with `P` open-loop unstable poles the closed loop has
#' \eqn{Z = N + P} unstable poles. The winding number is obtained from the
#' argument principle, \eqn{N = Z - P}, by counting the unstable roots of
#' the closed-loop and open-loop characteristic polynomials, so it does not
#' depend on how finely the plotted curve is sampled. It is `NA` when
#' \eqn{1 + L} vanishes on the contour (a closed-loop pole on the stability
#' boundary, or \eqn{1 + D = 0}). `ctrl_nichols()`
#' returns phase against
#' magnitude in dB.
#'
#' @param sys A model.
#' @param w Frequencies (rad per time unit).
#' @return Objects of class `ctrl_bode`, `ctrl_nyquist` or `ctrl_nichols`
#'   with a `data` data frame and a `plot()` method.
#' @examples
#' L <- ctrl_tf(10, c(1, 3, 3, 1))
#' b <- ctrl_bode(L)
#' head(b$data)
#' ny <- ctrl_nyquist(L)
#' ny$encirclements
#' @export
ctrl_bode <- function(sys, w = NULL) {
  fr <- ctrl_freqresp(sys, w)
  rows <- list()
  p <- dim(fr$response)[1]
  m <- dim(fr$response)[2]
  for (i in seq_len(p)) for (j in seq_len(m)) {
    g <- fr$response[i, j, ]
    ph <- unwrap_phase(Arg(g)) * 180 / pi
    rows[[length(rows) + 1L]] <- data.frame(
      w = fr$w, magnitude = Mod(g), magnitude_db = 20 * log10(Mod(g)),
      phase = ph, input = fr$input_names[j], output = fr$output_names[i],
      stringsAsFactors = FALSE)
  }
  data <- do.call(rbind, rows)
  structure(list(data = data, dt = fr$dt), class = "ctrl_bode")
}

#' @rdname ctrl_bode
#' @export
ctrl_nyquist <- function(sys, w = NULL) {
  fr <- ctrl_freqresp(sys, w)
  wpos <- fr$w
  rows <- list()
  p <- dim(fr$response)[1]
  m <- dim(fr$response)[2]
  for (i in seq_len(p)) for (j in seq_len(m)) {
    g <- fr$response[i, j, ]
    rows[[length(rows) + 1L]] <- data.frame(
      w = c(-rev(wpos), wpos), re = c(Re(rev(Conj(g))), Re(g)),
      im = c(Im(rev(Conj(g))), Im(g)), input = fr$input_names[j],
      output = fr$output_names[i], stringsAsFactors = FALSE)
  }
  data <- do.call(rbind, rows)
  enc <- NA_integer_
  if (p == 1L && m == 1L) {
    st <- ctrl_stability(sys)
    if (!length(st$boundary_poles)) {
      enc <- nyquist_encirclements(sys)
    }
  }
  structure(list(data = data, dt = fr$dt, encirclements = enc,
                 open_loop_unstable = if (p == 1L && m == 1L)
                   length(ctrl_stability(sys)$unstable_poles) else NA_integer_),
            class = "ctrl_nyquist")
}

# Clockwise encirclements of -1 by L along the Nyquist contour. By the
# argument principle N = Z - P, with Z and P the numbers of unstable roots
# of den + num (closed loop) and den (open loop) of L = num / den. Counting
# roots is exact up to their classification, whereas a sampled contour can
# miss a full turn between samples. NA when 1 + L vanishes on the contour.
nyquist_encirclements <- function(sys) {
  tf <- if (is_ctrl_tf(sys)) sys else if (inherits(sys, "ctrl_zpk")) ctrl_zpk2tf(sys) else
    ctrl_ss2tf(as_ss(sys), minimal = FALSE)
  num <- tf$num[[1]]
  den <- tf$den[[1]]
  L <- max(length(num), length(den))
  num <- c(rep(0, L - length(num)), num)
  den <- c(rep(0, L - length(den)), den)
  cl <- num + den
  # 1 + D = 0: the closed-loop polynomial loses degree, -1 is reached at
  # infinite (continuous) frequency
  if (abs(cl[1]) <= 1e-10 * max(abs(den[1]), abs(num[1]))) return(NA_integer_)
  unstable <- function(r, strict) {
    if (sys$dt > 0) {
      d <- Mod(r) - 1
      tol <- 1e-8
    } else {
      d <- Re(r)
      tol <- 1e-8 * pmax(1, Mod(r))
    }
    if (strict && any(abs(d) <= tol)) return(NA_integer_)
    sum(d > tol)
  }
  Z <- unstable(poly_roots(cl), TRUE)
  if (is.na(Z)) return(NA_integer_)
  P <- unstable(poly_roots(den), FALSE)
  as.integer(Z - P)
}

#' @rdname ctrl_bode
#' @export
ctrl_nichols <- function(sys, w = NULL) {
  b <- ctrl_bode(sys, w)
  structure(list(data = b$data, dt = b$dt), class = "ctrl_nichols")
}

siso_loop_tf <- function(L) {
  L <- as_system(L)
  if (is_ctrl_ss(L)) L <- ctrl_ss2tf(L)
  if (length(L$num) != 1L) {
    ctrl_abort("Classical gain and phase margins are defined for SISO loops ",
               "only; MIMO robustness needs disk margins or structured ",
               "singular values.")
  }
  L
}

# Polynomial in w (ascending coefficients) for Re and Im of p(jw), with p in
# descending powers of s.
jw_parts <- function(p) {
  a <- rev(p)
  k <- seq_along(a) - 1L
  re <- ifelse(k %% 2L == 0L, a * (-1)^(k %/% 2L), 0)
  im <- ifelse(k %% 2L == 1L, a * (-1)^((k - 1L) %/% 2L), 0)
  list(re = re, im = im)
}

asc_mul <- function(a, b) rev(poly_mul(rev(a), rev(b)))

asc_add <- function(a, b) {
  n <- max(length(a), length(b))
  c(a, rep(0, n - length(a))) + c(b, rep(0, n - length(b)))
}

positive_real_roots <- function(asc, tol = 1e-6) {
  asc <- asc[seq_len(max(which(abs(asc) > 0), 0))]
  if (length(asc) < 2L) return(numeric())
  r <- poly_roots(rev(asc))
  keep <- abs(Im(r)) <= tol * pmax(1, Mod(r)) & Re(r) > 0
  sort(unique(Re(r[keep])))
}

# Unit-circle roots of the Laurent polynomial given by coefficients of
# powers -d..d; returns angles in (0, pi].
circle_roots <- function(laurent, tol = 1e-5) {
  asc <- laurent
  nz <- which(abs(asc) > 0)
  if (!length(nz)) return(numeric())
  asc <- asc[min(nz):max(nz)]
  if (length(asc) < 2L) return(numeric())
  r <- poly_roots(rev(asc))
  keep <- abs(Mod(r) - 1) <= tol
  th <- abs(Arg(r[keep]))
  sort(unique(round(th[th > 1e-12], 12)))
}

laurent_prod <- function(a, b, d) {
  # a(z) * b(1/z) for ascending coefficient vectors a and b; return the
  # coefficients of powers -d..d
  out <- rep(0, 2L * d + 1L)
  for (i in seq_along(a)) for (j in seq_along(b)) {
    pw <- (i - 1L) - (j - 1L)
    out[pw + d + 1L] <- out[pw + d + 1L] + a[i] * b[j]
  }
  out
}

polish_root <- function(f, w, rel = 1e-3) {
  lo <- w * (1 - rel)
  hi <- w * (1 + rel)
  flo <- f(lo)
  fhi <- f(hi)
  if (is.finite(flo) && is.finite(fhi) && flo * fhi < 0) {
    return(stats::uniroot(f, c(lo, hi), tol = 1e-14 * max(1, w))$root)
  }
  w
}

#' Gain, phase and delay margins
#'
#' Classical stability margins of a SISO loop transfer function `L` (the
#' loop is closed with unity negative feedback). Gain crossover frequencies
#' (\eqn{|L| = 1}) and phase crossover frequencies (\eqn{L} real and
#' negative) are computed as roots of polynomials in \eqn{\omega} (or on the
#' unit circle for discrete models), then polished with a bracketing root
#' finder, so they do not depend on a frequency grid.
#'
#' The gain margin at a phase crossover \eqn{\omega_{pc}} is
#' \eqn{1/|L(j\omega_{pc})|}; the phase margin at a gain crossover
#' \eqn{\omega_{gc}} is \eqn{180 + \angle L(j\omega_{gc})} degrees,
#' wrapped to \eqn{(-180, 180]}; the delay margin is the smallest
#' \eqn{(\mathrm{PM}_i \bmod 360)/\omega_{gc,i}} over all gain crossovers
#' (radians over rad per time unit), the smallest added delay that brings a
#' crossover onto the critical point. The reported gain and phase margins
#' are the most critical ones (gain margin closest to 1 on a log scale,
#' smallest phase margin). `closed_loop_stable` is computed on the given
#' realization, so unstable hidden modes of a state-space loop are not
#' canceled away. Margins describe distance
#' to instability only if the nominal closed loop is stable, which is
#' reported in `closed_loop_stable`. Margins are SISO quantities; for MIMO
#' loops use loop-at-a-time disk margins.
#'
#' @param L A SISO loop transfer function (`ctrl_tf`, `ctrl_ss` or
#'   `ctrl_zpk`).
#' @return A list of class `ctrl_margin` with `gain_margin` (absolute),
#'   `gain_margin_db`, `phase_margin` (degrees), `delay_margin`,
#'   `phase_crossover` and `gain_crossover` frequencies, the tables `all_gm`
#'   and `all_pm`, and `closed_loop_stable`.
#' @examples
#' L <- ctrl_tf(10, c(1, 3, 3, 1))
#' ctrl_margin(L)
#' @export
ctrl_margin <- function(L) {
  L_orig <- as_system(L)
  L <- siso_loop_tf(L)
  num <- L$num[[1]]
  den <- L$den[[1]]
  dt <- L$dt
  evalL <- function(w) {
    z <- if (dt > 0) exp(1i * w * dt) else 1i * w
    poly_eval(num, z) / poly_eval(den, z)
  }
  if (dt == 0) {
    pn <- jw_parts(num)
    pd <- jw_parts(den)
    mag_poly <- asc_add(asc_add(asc_mul(pn$re, pn$re), asc_mul(pn$im, pn$im)),
                        -asc_add(asc_mul(pd$re, pd$re), asc_mul(pd$im, pd$im)))
    ph_poly <- asc_add(asc_mul(pn$im, pd$re), -asc_mul(pn$re, pd$im))
    wgc <- positive_real_roots(mag_poly)
    wpc <- positive_real_roots(ph_poly)
  } else {
    a <- rev(num)
    b <- rev(den)
    d <- max(length(a), length(b)) - 1L
    mag_l <- laurent_prod(a, a, d) - laurent_prod(b, b, d)
    ph_l <- laurent_prod(a, b, d) - laurent_prod(b, a, d)
    wgc <- circle_roots(mag_l) / dt
    wpc <- circle_roots(ph_l) / dt
    # the Nyquist frequency is always a phase crossover; drop near copies
    wpc <- c(wpc[abs(wpc * dt - pi) > 1e-9], pi / dt)
  }
  wgc <- vapply(wgc, function(w) polish_root(function(x) Mod(evalL(x)) - 1, w),
                0)
  wpc <- vapply(wpc, function(w) {
    if (dt > 0 && abs(w - pi / dt) < 1e-12) return(w)
    polish_root(function(x) Im(evalL(x)), w)
  }, 0)
  # include w = 0 when L(0) is finite, real and negative
  L0 <- suppressWarnings(evalL(0))
  if (is.finite(Mod(L0)) && Re(L0) < 0 && abs(Im(L0)) < 1e-12 * Mod(L0)) {
    wpc <- c(0, wpc)
  }
  gm_tab <- data.frame(frequency = numeric(), gain_margin = numeric())
  # a biproper continuous loop with L(j Inf) real and negative crosses the
  # negative real axis at infinite frequency
  if (dt == 0) {
    nn <- strip_leading_zeros(num)
    dd <- strip_leading_zeros(den)
    if (length(nn) == length(dd) && nn[1] / dd[1] < 0) {
      gm_tab <- data.frame(frequency = Inf, gain_margin = abs(dd[1] / nn[1]))
    }
  }
  for (w in unique(wpc)) {
    v <- evalL(w)
    if (!is.finite(Mod(v))) next
    if (Re(v) < 0 && abs(Im(v)) <= 1e-6 * Mod(v)) {
      gm_tab <- rbind(gm_tab, data.frame(frequency = w, gain_margin = 1 / Mod(v)))
    }
  }
  pm_tab <- data.frame(frequency = numeric(), phase_margin = numeric())
  for (w in unique(wgc)) {
    v <- evalL(w)
    ph <- Arg(v) * 180 / pi
    pm <- ((180 + ph + 180) %% 360) - 180
    if (pm <= -180) pm <- pm + 360
    pm_tab <- rbind(pm_tab, data.frame(frequency = w, phase_margin = pm))
  }
  gm <- Inf
  wpc_best <- NA_real_
  if (nrow(gm_tab)) {
    k <- which.min(abs(log(gm_tab$gain_margin)))
    gm <- gm_tab$gain_margin[k]
    wpc_best <- gm_tab$frequency[k]
  }
  pm <- Inf
  wgc_best <- NA_real_
  dm <- Inf
  if (nrow(pm_tab)) {
    k <- which.min(pm_tab$phase_margin)
    pm <- pm_tab$phase_margin[k]
    wgc_best <- pm_tab$frequency[k]
    # a delay adds phase lag, so the crossover at w_i is reached after a
    # lag of (PM_i mod 360) degrees; the delay margin is the smallest
    # such lag over all gain crossovers
    pm_tab$delay_margin <- ifelse(pm_tab$frequency > 0,
                                  ((pm_tab$phase_margin %% 360) * pi / 180) /
                                    pm_tab$frequency, 0)
    dm <- min(pm_tab$delay_margin)
  }
  # use the original realization: a minimal transfer function would hide
  # unstable uncontrollable or unobservable modes of the loop
  cl <- tryCatch(ctrl_is_stable(ctrl_feedback(L_orig, 1)), error = function(e) NA)
  structure(list(gain_margin = gm, gain_margin_db = 20 * log10(gm),
                 phase_margin = pm, delay_margin = dm,
                 phase_crossover = wpc_best, gain_crossover = wgc_best,
                 all_gm = gm_tab, all_pm = pm_tab, closed_loop_stable = cl,
                 dt = dt),
            class = "ctrl_margin")
}

#' @export
print.ctrl_margin <- function(x, digits = 4, ...) {
  cat("Stability margins (unity negative feedback)\n")
  cat("  gain margin: ", signif(x$gain_margin, digits), " (",
      signif(x$gain_margin_db, digits), " dB) at w = ",
      signif(x$phase_crossover, digits), "\n", sep = "")
  cat("  phase margin: ", signif(x$phase_margin, digits), " deg at w = ",
      signif(x$gain_crossover, digits), "\n", sep = "")
  cat("  delay margin: ", signif(x$delay_margin, digits), "\n", sep = "")
  cat("  nominal closed loop stable: ", x$closed_loop_stable, "\n", sep = "")
  invisible(x)
}

#' Bandwidth and resonant peak
#'
#' `ctrl_bandwidth()` returns the first frequency at which the gain of a
#' SISO model drops `drop_db` (default 3 dB) below its DC gain; it is `NA`
#' when the gain never drops that far on the analyzed range and an error
#' when the DC gain is infinite. `ctrl_peak()` returns the peak gain over
#' frequency, which equals the H-infinity norm for stable models.
#'
#' @param sys A SISO model.
#' @param drop_db Gain drop in dB.
#' @return `ctrl_bandwidth()`: a number. `ctrl_peak()`: a list with `gain`
#'   and `frequency`.
#' @examples
#' G <- ctrl_tf(1, c(1, 0.5, 1))
#' ctrl_bandwidth(G)
#' ctrl_peak(G)
#' @export
ctrl_bandwidth <- function(sys, drop_db = 3) {
  sys <- siso_loop_tf(sys)
  g0 <- abs(ctrl_dcgain(sys)[1, 1])
  target <- g0 * 10^(-drop_db / 20)
  w <- ctrl_freq_grid(sys, 2000L)
  mag <- Mod(eval_system(sys, if (sys$dt > 0) exp(1i * w * sys$dt) else
    1i * w)[1, 1, ])
  k <- which(mag < target)[1]
  if (is.na(k)) return(NA_real_)
  if (k == 1L) return(w[1])
  f <- function(x) {
    z <- if (sys$dt > 0) exp(1i * x * sys$dt) else 1i * x
    Mod(eval_system(sys, z)[1, 1, ]) - target
  }
  stats::uniroot(f, c(w[k - 1L], w[k]), tol = 1e-12)$root
}

#' @rdname ctrl_bandwidth
#' @export
ctrl_peak <- function(sys) {
  n <- ctrl_norm(sys, "inf")
  list(gain = as.numeric(n), frequency = attr(n, "frequency"))
}
