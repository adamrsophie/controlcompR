#' Plot methods
#'
#' All plots are built with 'ggplot2' and returned as `ggplot` objects, so
#' they can be modified or saved and work in noninteractive sessions.
#' Discrete time signals are drawn as steps (sample and hold); continuous
#' signals as lines. Axis labels state the time domain and units
#' (rad per time unit for frequencies).
#'
#' @param x An object returned by a `ctrl_*` analysis function.
#' @param what For responses: `"outputs"`, `"states"`, `"inputs"` or
#'   `"all"`.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @name ctrl_plots
#' @examples
#' G <- ctrl_tf(4, c(1, 1.2, 4))
#' plot(ctrl_step(G))
#' plot(ctrl_bode(G))
#' plot(ctrl_nyquist(G))
#' plot(ctrl_pzmap(G))
#' plot(ctrl_rlocus(ctrl_tf(1, c(1, 3, 2, 0))))
NULL

time_label <- function(dt) {
  if (dt > 0) paste0("time (sampled, dt = ", format(dt), ")") else "time"
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_response <- function(x, what = c("outputs", "states", "inputs", "all"),
                               ...) {
  what <- match.arg(what)
  df <- as.data.frame(x)
  kinds <- switch(what, outputs = "output", states = "state", inputs = "input",
                  all = c("output", "state", "input"))
  df <- df[df$kind %in% kinds, , drop = FALSE]
  df$panel <- paste(df$kind, df$signal, sep = ": ")
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value,
                                        color = .data$case))
  g <- g + if (x$dt > 0) ggplot2::geom_step() else ggplot2::geom_line()
  g + ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::labs(x = time_label(x$dt), y = NULL, color = NULL,
                  title = paste(tools_title(x$type), "response")) +
    ggplot2::theme_bw()
}

tools_title <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))

#' @rdname ctrl_plots
#' @export
plot.ctrl_bode <- function(x, ...) {
  df <- x$data
  df$channel <- paste(df$input, "->", df$output)
  long <- rbind(
    data.frame(w = df$w, value = df$magnitude_db, panel = "magnitude (dB)",
               channel = df$channel),
    data.frame(w = df$w, value = df$phase, panel = "phase (deg)",
               channel = df$channel)
  )
  long$panel <- factor(long$panel, c("magnitude (dB)", "phase (deg)"))
  g <- ggplot2::ggplot(long, ggplot2::aes(x = .data$w, y = .data$value,
                                          color = .data$channel)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_log10() +
    ggplot2::facet_grid(panel ~ ., scales = "free_y") +
    ggplot2::labs(x = "frequency (rad per time unit)", y = NULL, color = NULL,
                  title = "Bode diagram") +
    ggplot2::theme_bw()
  if (x$dt > 0) {
    g <- g + ggplot2::geom_vline(xintercept = pi / x$dt, linetype = "dashed")
  }
  g
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_nyquist <- function(x, ...) {
  df <- x$data
  df$branch <- ifelse(df$w >= 0, "w > 0", "w < 0")
  df$channel <- paste(df$input, "->", df$output)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$re, y = .data$im,
                                   linetype = .data$branch,
                                   color = .data$channel)) +
    ggplot2::geom_path() +
    ggplot2::annotate("point", x = -1, y = 0, shape = 3, size = 3,
                      color = "red") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "real", y = "imaginary", linetype = NULL, color = NULL,
                  title = "Nyquist diagram",
                  subtitle = if (!is.na(x$encirclements))
                    paste("clockwise encirclements of -1:", x$encirclements)) +
    ggplot2::theme_bw()
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_nichols <- function(x, ...) {
  df <- x$data
  df$channel <- paste(df$input, "->", df$output)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$phase, y = .data$magnitude_db,
                                   color = .data$channel)) +
    ggplot2::geom_path() +
    ggplot2::annotate("point", x = -180, y = 0, shape = 3, size = 3,
                      color = "red") +
    ggplot2::labs(x = "phase (deg)", y = "magnitude (dB)", color = NULL,
                  title = "Nichols chart") +
    ggplot2::theme_bw()
}

#' Pole-zero map data
#'
#' @param sys A model.
#' @return An object of class `ctrl_pzmap` with `poles`, `zeros` and `dt`.
#' @examples
#' ctrl_pzmap(ctrl_tf(c(1, 2), c(1, 2, 5)))
#' @export
ctrl_pzmap <- function(sys) {
  z <- tryCatch(ctrl_zeros(sys), error = function(e) complex())
  structure(list(poles = ctrl_poles(sys), zeros = z, dt = sys$dt),
            class = "ctrl_pzmap")
}

#' @export
print.ctrl_pzmap <- function(x, digits = 4, ...) {
  cat("Poles: ", paste(format(signif(x$poles, digits)), collapse = ", "),
      "\nZeros: ", if (length(x$zeros)) paste(format(signif(x$zeros, digits)),
                                              collapse = ", ") else "none",
      "\n", sep = "")
  invisible(x)
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_pzmap <- function(x, ...) {
  df <- rbind(
    data.frame(re = Re(x$poles), im = Im(x$poles), type = rep("pole", length(x$poles))),
    data.frame(re = Re(x$zeros), im = Im(x$zeros), type = rep("zero", length(x$zeros)))
  )
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$re, y = .data$im,
                                        shape = .data$type)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, color = "grey60") +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_shape_manual(values = c(pole = 4, zero = 1)) +
    ggplot2::labs(x = "real", y = "imaginary", shape = NULL,
                  title = "Pole-zero map") +
    ggplot2::theme_bw()
  if (x$dt > 0) {
    th <- seq(0, 2 * pi, length.out = 361)
    g <- g + ggplot2::annotate("path", x = cos(th), y = sin(th),
                               linetype = "dashed", color = "grey40") +
      ggplot2::coord_equal()
  }
  g
}

#' Root locus
#'
#' Closed-loop poles of \eqn{1 + k L(s) = 0} for a SISO loop as the gain
#' `k` varies, computed as the roots of \eqn{d(s) + k\,n(s)}. Branches are
#' tracked by matching each root to the nearest root at the previous gain.
#'
#' @param L SISO loop transfer function.
#' @param gains Optional nonnegative gains; by default a logarithmic grid
#'   chosen so that the poles approach the open-loop zeros or asymptotes.
#' @return An object of class `ctrl_rlocus` with `gains` and a complex
#'   matrix `roots` (gains x branches).
#' @examples
#' rl <- ctrl_rlocus(ctrl_tf(1, c(1, 3, 2, 0)))
#' dim(rl$roots)
#' @export
ctrl_rlocus <- function(L, gains = NULL) {
  L <- siso_loop_tf(L)
  num <- L$num[[1]]
  den <- L$den[[1]]
  n <- length(den) - 1L
  if (n < 1L) ctrl_abort("The loop has no poles.")
  if (is.null(gains)) {
    scale <- max(abs(den)) / max(abs(num))
    gains <- c(0, scale * 10^seq(-4, 4, length.out = 400))
  }
  if (any(gains < 0)) ctrl_abort("`gains` must be nonnegative.")
  gains <- sort(gains)
  R <- matrix(NA_complex_, length(gains), n)
  prev <- NULL
  for (k in seq_along(gains)) {
    cp <- poly_add(den, gains[k] * num)
    r <- poly_roots(cp)
    r <- c(r, rep(NA_complex_, n - length(r)))[seq_len(n)]
    if (!is.null(prev)) {
      ord <- integer(0)
      avail <- seq_len(n)
      for (j in seq_len(n)) {
        if (is.na(prev[j])) {
          ord <- c(ord, avail[1])
          avail <- avail[-1]
          next
        }
        d <- Mod(r[avail] - prev[j])
        d[is.na(d)] <- Inf
        pick <- avail[which.min(d)]
        ord <- c(ord, pick)
        avail <- setdiff(avail, pick)
      }
      r <- r[ord]
    }
    R[k, ] <- r
    prev <- r
  }
  structure(list(gains = gains, roots = R, poles = poly_roots(den),
                 zeros = poly_roots(num), dt = L$dt),
            class = "ctrl_rlocus")
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_rlocus <- function(x, ...) {
  df <- do.call(rbind, lapply(seq_len(ncol(x$roots)), function(j) {
    data.frame(re = Re(x$roots[, j]), im = Im(x$roots[, j]),
               branch = factor(j), gain = x$gains)
  }))
  df <- df[is.finite(df$re), , drop = FALSE]
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$re, y = .data$im,
                                        color = .data$branch)) +
    ggplot2::geom_path() +
    ggplot2::annotate("point", x = Re(x$poles), y = Im(x$poles), shape = 4,
                      size = 3) +
    ggplot2::labs(x = "real", y = "imaginary", color = "branch",
                  title = "Root locus") +
    ggplot2::theme_bw()
  if (length(x$zeros)) {
    g <- g + ggplot2::annotate("point", x = Re(x$zeros), y = Im(x$zeros),
                               shape = 1, size = 3)
  }
  if (x$dt > 0) {
    th <- seq(0, 2 * pi, length.out = 361)
    g <- g + ggplot2::annotate("path", x = cos(th), y = sin(th),
                               linetype = "dashed", color = "grey40")
  }
  g
}
