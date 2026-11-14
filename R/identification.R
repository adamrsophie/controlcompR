#' Identification data
#'
#' Store sampled input-output data for system identification.
#'
#' @param y Output data (vector or `T x p` matrix).
#' @param u Input data (vector or `T x m` matrix).
#' @param dt Sample interval.
#' @param detrend Remove the sample means of `y` and `u` (stored in
#'   `offsets` so predictions can be mapped back), or a list with `y` and `u`
#'   offsets to subtract, typically `train$offsets`. Models estimated on
#'   detrended data describe deviations from the offsets, so validation data
#'   must be shifted by the training offsets, not by their own means.
#' @return An object of class `ctrl_iddata`.
#' @examples
#' d <- ctrl_iddata(y = sin(1:50 / 5), u = cos(1:50 / 5), dt = 0.1)
#' d
#' @export
ctrl_iddata <- function(y, u, dt = 1, detrend = FALSE) {
  Y <- if (is.null(dim(y))) matrix(as.numeric(y), ncol = 1) else as.matrix(y)
  U <- if (is.null(dim(u))) matrix(as.numeric(u), ncol = 1) else as.matrix(u)
  storage.mode(Y) <- "double"
  storage.mode(U) <- "double"
  if (nrow(Y) != nrow(U)) ctrl_abort("`y` and `u` must have the same number of samples.")
  if (anyNA(Y) || anyNA(U)) ctrl_abort("Identification data must not contain NA.")
  check_scalar_number(dt, "dt", positive = TRUE)
  off <- list(y = rep(0, ncol(Y)), u = rep(0, ncol(U)))
  if (is.list(detrend)) {
    oy <- as.numeric(detrend$y)
    ou <- as.numeric(detrend$u)
    if (length(oy) != ncol(Y) || length(ou) != ncol(U)) {
      ctrl_abort("`detrend` offsets must have one value per output and input.")
    }
    off <- list(y = oy, u = ou)
    Y <- sweep(Y, 2, off$y)
    U <- sweep(U, 2, off$u)
  } else if (isTRUE(detrend)) {
    off <- list(y = colMeans(Y), u = colMeans(U))
    Y <- sweep(Y, 2, off$y)
    U <- sweep(U, 2, off$u)
  }
  colnames(Y) <- colnames(Y) %||% default_names("y", ncol(Y))
  colnames(U) <- colnames(U) %||% default_names("u", ncol(U))
  structure(list(y = Y, u = U, dt = dt, offsets = off), class = "ctrl_iddata")
}

#' @export
print.ctrl_iddata <- function(x, ...) {
  cat("Identification data: ", nrow(x$y), " samples, ", ncol(x$u), " input(s), ",
      ncol(x$y), " output(s), dt = ", x$dt, "\n", sep = "")
  invisible(x)
}

#' Split identification data into estimation and validation parts
#'
#' @param data A `ctrl_iddata` object.
#' @param fraction Fraction of samples used for estimation (the first
#'   part; the remainder is kept for validation, preserving time order).
#' @return A list with `estimation` and `validation` data sets.
#' @examples
#' d <- ctrl_iddata(rnorm(100), rnorm(100))
#' parts <- ctrl_idsplit(d, 0.7)
#' @export
ctrl_idsplit <- function(data, fraction = 0.7) {
  check_iddata(data)
  Tn <- nrow(data$y)
  k <- floor(fraction * Tn)
  if (k < 2L || k > Tn - 2L) ctrl_abort("`fraction` leaves too few samples in a part.")
  sub <- function(idx) {
    out <- data
    out$y <- data$y[idx, , drop = FALSE]
    out$u <- data$u[idx, , drop = FALSE]
    out
  }
  list(estimation = sub(seq_len(k)), validation = sub((k + 1L):Tn))
}

check_iddata <- function(data) {
  if (!inherits(data, "ctrl_iddata")) ctrl_abort("`data` must be a `ctrl_iddata` object.")
  invisible(data)
}

#' ARX model estimation
#'
#' Estimate the single-output ARX model
#' \deqn{y_k + a_1 y_{k-1} + \cdots + a_{n_a} y_{k-n_a} =
#' \sum_{i=1}^{m} \sum_{j=0}^{n_b - 1} b_{ij} u_{i, k - n_k - j} + e_k}
#' by linear least squares (QR). The convention for `nb` and `nk` matches
#' the `sysid` package: `nb` coefficients starting at delay `nk`. Standard
#' errors assume white equation errors; they are not valid when the true
#' noise is colored (use an output-error or Box-Jenkins model from
#' `sysid` through [ctrl_from_sysid()] in that case).
#'
#' @param data A `ctrl_iddata` object with one output.
#' @param na Number of output (denominator) coefficients.
#' @param nb Number of input coefficients per input.
#' @param nk Input delay in samples.
#' @return An object of class `c("ctrl_arx", "ctrl_idmodel")` with `sys`
#'   (a discrete `ctrl_tf` from inputs to the output), `a`, `b`, `vcov`,
#'   `noise_var`, `residuals` and information criteria.
#' @examples
#' set.seed(1)
#' u <- rnorm(300)
#' y <- as.numeric(stats::filter(0.5 * c(0, u[-300]), 0.8, "recursive")) +
#'   rnorm(300, sd = 0.05)
#' m <- ctrl_arx(ctrl_iddata(y, u), na = 1, nb = 1, nk = 1)
#' m$a
#' m$b
#' @export
ctrl_arx <- function(data, na, nb, nk = 1) {
  check_iddata(data)
  if (ncol(data$y) != 1L) ctrl_abort("`ctrl_arx()` estimates single-output models; use ctrl_n4sid() for MIMO data.")
  for (nm in c("na", "nb", "nk")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || v < 0 || v != round(v)) {
      ctrl_abort("`", nm, "` must be a nonnegative integer.")
    }
  }
  if (nb < 1L) ctrl_abort("`nb` must be at least 1.")
  y <- data$y[, 1]
  U <- data$u
  m <- ncol(U)
  Tn <- length(y)
  k0 <- max(na, nk + nb - 1L) + 1L
  rows <- k0:Tn
  if (length(rows) <= na + m * nb) ctrl_abort("Too few samples for the requested orders.")
  X <- cbind(
    if (na) vapply(seq_len(na), function(i) -y[rows - i], numeric(length(rows))),
    do.call(cbind, lapply(seq_len(m), function(i) {
      vapply(seq_len(nb), function(j) U[rows - nk - j + 1L, i], numeric(length(rows)))
    }))
  )
  X <- matrix(X, nrow = length(rows))
  q <- qr(X)
  if (q$rank < ncol(X)) {
    ctrl_abort("The ARX regression is rank deficient; the input is not ",
               "persistently exciting for these orders.")
  }
  theta <- qr.coef(q, y[rows])
  res <- y[rows] - drop(X %*% theta)
  npar <- length(theta)
  N <- length(rows)
  sigma2 <- sum(res^2) / (N - npar)
  R <- qr.R(q)
  Rinv <- backsolve(R, diag(npar))
  vc <- sigma2 * Rinv %*% t(Rinv)
  vc <- vc[order(q$pivot), order(q$pivot), drop = FALSE]
  a <- theta[seq_len(na)]
  b <- matrix(theta[na + seq_len(m * nb)], nb, m)
  den <- c(1, a)
  num <- matrix(list(), 1, m)
  dn <- matrix(list(), 1, m)
  for (i in seq_len(m)) {
    nq <- c(rep(0, nk), b[, i])
    L <- max(length(nq), length(den))
    num[[1, i]] <- c(nq, rep(0, L - length(nq)))
    dn[[1, i]] <- c(den, rep(0, L - length(den)))
  }
  sys <- ctrl_tf(num, dn, dt = data$dt, input_names = colnames(data$u),
                 output_names = colnames(data$y))
  ll <- -0.5 * N * (log(2 * pi * sum(res^2) / N) + 1)
  structure(list(sys = sys, a = a, b = b, na = na, nb = nb, nk = nk,
                 vcov = vc, noise_var = sigma2, residuals = res,
                 loglik = ll, aic = -2 * ll + 2 * npar,
                 bic = -2 * ll + log(N) * npar, npar = npar, nobs = N,
                 method = "arx", dt = data$dt),
            class = c("ctrl_arx", "ctrl_idmodel"))
}

#' Subspace state-space identification (PO-MOESP)
#'
#' Estimate a discrete time state-space model
#' \eqn{x_{k+1} = A x_k + B u_k}, \eqn{y_k = C x_k + D u_k + v_k} from
#' input-output data with the past-output MOESP subspace method (Verhaegen,
#' 1994). Future outputs are projected on future inputs and instrumented by
#' past inputs and outputs (an LQ factorization of stacked block Hankel
#' matrices); the singular values of the projection reveal the order, the
#' extended observability matrix gives `A` and `C` (shift invariance), and
#' `B`, `D` and the initial state are then estimated by linear least
#' squares on the simulated output.
#'
#' The instruments make the estimates consistent for white or colored
#' output noise in open loop. Closed-loop data (inputs correlated with past
#' noise) bias the estimates; this is not corrected. The model has no
#' noise model (output-error form).
#'
#' @param data A `ctrl_iddata` object.
#' @param order Model order; `NULL` chooses the order at the largest gap of
#'   the log singular values.
#' @param horizon Number of block rows `i` (must exceed the order).
#' @param feedthrough Estimate `D` (default `FALSE`, so `D = 0`).
#' @return An object of class `c("ctrl_n4sid", "ctrl_idmodel")` with `sys`
#'   (a `ctrl_ss`), `x0`, `singular_values`, `order` and `residuals`.
#' @references Verhaegen, M. (1994). Identification of the deterministic
#'   part of MIMO state space models given in innovations form from
#'   input-output data. *Automatica*, 30(1), 61-74.
#' @examples
#' set.seed(2)
#' true <- ctrl_ss(matrix(c(0.7, 0.2, -0.2, 0.7), 2), c(1, 0), c(1, 1), dt = 1)
#' u <- rnorm(400)
#' y <- ctrl_simulate(true, input = u, times = 0:399)$y[, 1, 1] + rnorm(400, sd = 0.02)
#' m <- ctrl_n4sid(ctrl_iddata(y, u), order = 2)
#' ctrl_poles(m$sys)
#' @export
ctrl_n4sid <- function(data, order = NULL, horizon = NULL, feedthrough = FALSE) {
  check_iddata(data)
  Y <- data$y
  U <- data$u
  Tn <- nrow(Y)
  p <- ncol(Y)
  m <- ncol(U)
  i <- if (is.null(horizon)) max(2L, min(15L, floor(Tn / (4 * (m + p + 1))))) else
    as_count(horizon, "horizon")
  if (!is.null(order) && order >= i * p) {
    ctrl_abort("`horizon` must be larger than order / p; increase `horizon`.")
  }
  N <- Tn - 2L * i + 1L
  if (N < 2L * i * (m + p)) ctrl_abort("Too few samples for horizon ", i, ".")
  hank <- function(X, start) {
    k <- ncol(X)
    H <- matrix(0, i * k, N)
    for (r in seq_len(i)) H[(r - 1L) * k + seq_len(k), ] <- t(X[start + r - 1L + seq_len(N) - 1L, , drop = FALSE])
    H
  }
  Up <- hank(U, 1L)
  Uf <- hank(U, i + 1L)
  Yp <- hank(Y, 1L)
  Yf <- hank(Y, i + 1L)
  # the input must be persistently exciting of order 2i
  ur <- qr(t(rbind(Uf, Up)))$rank
  if (ur < 2L * i * m) {
    ctrl_abort("The input is not persistently exciting of order ", 2L * i,
               " (rank ", ur, " of ", 2L * i * m, "); the system cannot be ",
               "identified from these data. Use a richer input or a shorter horizon.")
  }
  M <- rbind(Uf, Up, Yp, Yf) / sqrt(N)
  # tol = 0 disables column pivoting, which would scramble the block structure
  qm <- qr(t(M), tol = 0)
  if (!identical(qm$pivot, seq_len(nrow(M)))) {
    ctrl_abort("Internal error: pivoted QR in subspace identification.")
  }
  L <- t(qr.R(qm))
  nu <- i * m
  nz <- i * (m + p)
  iy <- nu + nz + seq_len(i * p)
  L32 <- L[iy, nu + seq_len(nz), drop = FALSE]
  sv <- svd(L32)
  sval <- sv$d
  if (is.null(order)) {
    ls <- log(pmax(sval[seq_len(min(length(sval), i * p - 1L))], .Machine$double.eps))
    # with a single admissible order (i p = 2) there is no gap to compare
    order <- if (length(ls) < 2L) 1L else which.max(-diff(ls))
  }
  n <- as_count(order, "order")
  Gam <- sv$u[, seq_len(n), drop = FALSE] %*% diag(sqrt(sval[seq_len(n)]), n)
  C <- Gam[seq_len(p), , drop = FALSE]
  up <- Gam[seq_len((i - 1L) * p), , drop = FALSE]
  dn <- Gam[p + seq_len((i - 1L) * p), , drop = FALSE]
  A <- qr.solve(up, dn)
  # B, D and x0 by least squares on the simulated output
  bd <- ss_input_ls(A, C, U, Y, feedthrough)
  sys <- ctrl_ss(A, bd$B, C, bd$D, dt = data$dt, input_names = colnames(U),
                 output_names = colnames(Y))
  structure(list(sys = sys, x0 = bd$x0, singular_values = sval, order = n,
                 horizon = i, residuals = bd$residuals, method = "n4sid",
                 dt = data$dt, npar = n * (n + m + p) + if (feedthrough) p * m else 0,
                 nobs = Tn),
            class = c("ctrl_n4sid", "ctrl_idmodel"))
}

ss_input_ls <- function(A, C, U, Y, feedthrough) {
  n <- nrow(A)
  p <- nrow(C)
  m <- ncol(U)
  Tn <- nrow(Y)
  nb <- n * m
  nd <- if (feedthrough) p * m else 0L
  Phi <- matrix(0, Tn * p, n + nb + nd)
  Ok <- diag(n)
  Sk <- matrix(0, n, nb)
  for (k in seq_len(Tn)) {
    rows <- (k - 1L) * p + seq_len(p)
    Phi[rows, seq_len(n)] <- C %*% Ok
    Phi[rows, n + seq_len(nb)] <- C %*% Sk
    if (feedthrough) Phi[rows, n + nb + seq_len(nd)] <- kronecker(t(U[k, ]), diag(p))
    Sk <- A %*% Sk + kronecker(t(U[k, ]), diag(n))
    Ok <- A %*% Ok
  }
  yv <- as.vector(t(Y))
  th <- qr.coef(qr(Phi), yv)
  th[is.na(th)] <- 0
  res <- matrix(yv - drop(Phi %*% th), Tn, p, byrow = TRUE)
  list(x0 = th[seq_len(n)], B = matrix(th[n + seq_len(nb)], n, m),
       D = if (feedthrough) matrix(th[n + nb + seq_len(nd)], p, m) else matrix(0, p, m),
       residuals = res)
}

#' Convert a `sysid` model
#'
#' Convert a model estimated by the `sysid` package (an `estpoly` fit or
#' an `idpoly` object) to `controlcompR` transfer functions. The `sysid`
#' polynomial model
#' \deqn{A(q) y = \frac{B(q)}{F(q)} q^{-n_k} u + \frac{C(q)}{D(q)} e}
#' (polynomials in \eqn{q^{-1}}) gives the plant
#' \eqn{G(z) = z^{-n_k} B(z^{-1}) / (A(z^{-1}) F(z^{-1}))} and the noise
#' model \eqn{H(z) = C(z^{-1}) / (A(z^{-1}) D(z^{-1}))}.
#'
#' @param model An `estpoly` or `idpoly` object from `sysid`.
#' @return A list with `plant` and `noise` (`ctrl_tf`), `noise_var` and
#'   `type`.
#' @examples
#' if (requireNamespace("sysid", quietly = TRUE)) {
#'   set.seed(1)
#'   u <- rnorm(200)
#'   y <- as.numeric(stats::filter(0.5 * c(0, u[-200]), 0.8, "recursive")) +
#'     rnorm(200, sd = 0.05)
#'   fit <- sysid::arx(sysid::idframe(output = y, input = u, Ts = 1), c(1, 1, 1))
#'   ctrl_from_sysid(fit)$plant
#' }
#' @export
ctrl_from_sysid <- function(model) {
  sys <- if (inherits(model, "estpoly")) model$sys else model
  if (!inherits(sys, "idpoly")) ctrl_abort("`model` must be a sysid `estpoly` or `idpoly` object.")
  qinv_tf <- function(num, den, dt) {
    L <- max(length(num), length(den))
    ctrl_tf(c(num, rep(0, L - length(num))), c(den, rep(0, L - length(den))), dt = dt)
  }
  dt <- sys$Ts %||% 1
  B <- c(rep(0, sys$ioDelay %||% 0), sys$B)
  AF <- poly_mul(sys$A, sys$F1 %||% 1)
  AD <- poly_mul(sys$A, sys$D %||% 1)
  list(plant = qinv_tf(B, AF, dt), noise = qinv_tf(sys$C %||% 1, AD, dt),
       noise_var = sys$noiseVar, type = sys$type)
}

# Simulated and one-step predicted outputs of an identified model.
idmodel_outputs <- function(model, data) {
  Tn <- nrow(data$y)
  times <- (seq_len(Tn) - 1) * data$dt
  if (inherits(model, "ctrl_n4sid")) {
    # estimate only the initial state on these data, keeping the model's
    # A, B, C and D fixed
    s0 <- ctrl_simulate(model$sys, input = data$u, times = times)$y[, , 1]
    s0 <- matrix(s0, Tn)
    nx <- nrow(model$sys$A)
    O <- matrix(0, Tn * ncol(data$y), nx)
    Ak <- diag(nx)
    for (k in seq_len(Tn)) {
      O[(k - 1L) * ncol(data$y) + seq_len(ncol(data$y)), ] <- model$sys$C %*% Ak
      Ak <- model$sys$A %*% Ak
    }
    x0 <- qr.solve(O, as.numeric(t(data$y - s0)))
    sim <- ctrl_simulate(model$sys, input = data$u, times = times, x0 = x0)$y[, , 1]
    return(list(simulated = matrix(sim, Tn), predicted = NULL))
  }
  sys <- if (inherits(model, "ctrl_idmodel")) model$sys else as_system(model)
  sim <- ctrl_simulate(sys, input = data$u, times = times)$y[, , 1]
  pred <- NULL
  if (inherits(model, "ctrl_arx")) {
    y <- data$y[, 1]
    U <- data$u
    pred <- rep(NA_real_, Tn)
    for (k in seq_len(Tn)) {
      v <- 0
      for (i in seq_len(model$na)) v <- v - model$a[i] * if (k - i >= 1) y[k - i] else 0
      for (j in seq_len(ncol(U))) for (l in seq_len(model$nb)) {
        kk <- k - model$nk - l + 1L
        if (kk >= 1) v <- v + model$b[l, j] * U[kk, j]
      }
      pred[k] <- v
    }
    pred <- matrix(pred, Tn)
  }
  list(simulated = matrix(sim, Tn), predicted = pred)
}

nrmse_fit <- function(y, yhat) {
  vapply(seq_len(ncol(y)), function(j) {
    100 * (1 - sqrt(sum((y[, j] - yhat[, j])^2)) /
             sqrt(sum((y[, j] - mean(y[, j]))^2)))
  }, 0)
}

#' Compare identified models with data
#'
#' Compute normalized fit percentages,
#' \eqn{100 (1 - \|y - \hat y\| / \|y - \bar y\|)}, for the simulated
#' output (input only, the demanding test relevant for control design) and,
#' for ARX models, the one-step-ahead prediction. Use validation data that
#' were not used for estimation. A high prediction fit does not establish
#' that a model is adequate for feedback design.
#'
#' @param data A `ctrl_iddata` object (preferably validation data).
#' @param ... Identified models (`ctrl_idmodel` objects or plain models),
#'   optionally named.
#' @return A data frame with one row per model and output, with class
#'   `ctrl_idcompare` (the simulated outputs are attached as attribute
#'   `outputs` for plotting).
#' @examples
#' set.seed(1)
#' u <- rnorm(400)
#' y <- as.numeric(stats::filter(0.5 * c(0, u[-400]), 0.8, "recursive")) +
#'   rnorm(400, sd = 0.05)
#' parts <- ctrl_idsplit(ctrl_iddata(y, u))
#' m1 <- ctrl_arx(parts$estimation, 1, 1, 1)
#' m2 <- ctrl_arx(parts$estimation, 2, 2, 1)
#' ctrl_compare(parts$validation, first = m1, second = m2)
#' @export
ctrl_compare <- function(data, ...) {
  check_iddata(data)
  models <- list(...)
  nms <- names(models) %||% rep("", length(models))
  nms[nms == ""] <- paste0("model", which(nms == ""))
  rows <- list()
  outs <- list()
  for (k in seq_along(models)) {
    o <- idmodel_outputs(models[[k]], data)
    fs <- nrmse_fit(data$y, o$simulated)
    fp <- if (is.null(o$predicted)) rep(NA_real_, ncol(data$y)) else
      nrmse_fit(data$y, o$predicted)
    rows[[k]] <- data.frame(model = nms[k], output = colnames(data$y),
                            simulation_fit = fs, prediction_fit = fp)
    outs[[nms[k]]] <- o$simulated
  }
  out <- do.call(rbind, rows)
  attr(out, "outputs") <- outs
  attr(out, "data") <- data
  class(out) <- c("ctrl_idcompare", "data.frame")
  out
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_idcompare <- function(x, ...) {
  data <- attr(x, "data")
  outs <- attr(x, "outputs")
  t <- (seq_len(nrow(data$y)) - 1) * data$dt
  rows <- list(data.frame(time = t, value = data$y[, 1], series = "measured"))
  for (nm in names(outs)) {
    rows[[length(rows) + 1L]] <- data.frame(time = t, value = outs[[nm]][, 1],
                                            series = nm)
  }
  ggplot2::ggplot(do.call(rbind, rows),
                  ggplot2::aes(x = .data$time, y = .data$value,
                               color = .data$series)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "time", y = colnames(data$y)[1], color = NULL,
                  title = "Measured and simulated outputs") +
    ggplot2::theme_bw()
}

#' Residual analysis for identified models
#'
#' Compute the residuals of a model on data (one-step prediction errors for
#' ARX models, output errors otherwise) and test them: the autocorrelation
#' should be within the 99 percent band \eqn{\pm 2.58/\sqrt{N}} for white
#' residuals (Ljung-Box test reported), and the cross-correlation with the
#' input should vanish at nonnegative lags (residuals independent of current
#' and past inputs; `ccf` reports \eqn{\mathrm{cor}(e_{t+k}, u_t)}). All
#' outputs and input-output pairs are tested. Output-error
#' models are not expected to have white residuals; the cross-correlation
#' test is the relevant one for them.
#'
#' @param model An identified model.
#' @param data Data (preferably validation data).
#' @param lags Maximum lag.
#' @return A list of class `ctrl_resid` with `residuals`, `acf`,
#'   `ccf` (data frames with output, input, lag, value and band),
#'   `ljung_box` (p-value per output) and
#'   flags `white` and `independent`.
#' @examples
#' set.seed(1)
#' u <- rnorm(300)
#' y <- as.numeric(stats::filter(0.5 * c(0, u[-300]), 0.8, "recursive")) +
#'   rnorm(300, sd = 0.05)
#' d <- ctrl_iddata(y, u)
#' r <- ctrl_resid(ctrl_arx(d, 1, 1, 1), d)
#' r$white
#' @export
ctrl_resid <- function(model, data, lags = 20L) {
  check_iddata(data)
  o <- idmodel_outputs(model, data)
  yhat <- if (!is.null(o$predicted)) o$predicted else o$simulated
  E <- as.matrix(data$y) - as.matrix(yhat)
  U <- as.matrix(data$u)
  skip <- if (inherits(model, "ctrl_arx")) seq_len(max(model$na, model$nk + model$nb - 1L)) else integer()
  if (length(skip)) {
    E <- E[-skip, , drop = FALSE]
    U <- U[-skip, , drop = FALSE]
  }
  N <- nrow(E)
  band <- 2.58 / sqrt(N)
  fitdf <- if (inherits(model, "ctrl_arx")) model$na else 0
  acd <- ccd <- list()
  lb <- numeric(ncol(E))
  for (i in seq_len(ncol(E))) {
    e <- E[, i]
    ac <- stats::acf(e, lag.max = lags, plot = FALSE)$acf[-1]
    acd[[i]] <- data.frame(output = i, lag = seq_len(lags), value = ac,
                           band = band)
    lb[i] <- stats::Box.test(e, lag = lags, type = "Ljung-Box",
                             fitdf = fitdf)$p.value
    for (j in seq_len(ncol(U))) {
      # ccf(e, u) at lag k estimates cor(e[t + k], u[t]): k >= 0 measures the
      # dependence of the residual on current and past inputs
      cc <- stats::ccf(e, U[, j], lag.max = lags, plot = FALSE)
      ccd[[length(ccd) + 1L]] <- data.frame(output = i, input = j,
                                            lag = drop(cc$lag),
                                            value = drop(cc$acf), band = band)
    }
  }
  acd <- do.call(rbind, acd)
  ccd <- do.call(rbind, ccd)
  structure(list(residuals = if (ncol(E) == 1L) E[, 1] else E,
                 acf = acd, ccf = ccd, ljung_box = lb,
                 white = all(abs(acd$value) < band),
                 independent = all(abs(ccd$value[ccd$lag >= 0]) < band)),
            class = "ctrl_resid")
}

#' @export
print.ctrl_resid <- function(x, ...) {
  cat("Residual analysis (", NROW(x$residuals), " samples)\n", sep = "")
  cat("  autocorrelation within 99% band: ", x$white, " (Ljung-Box p = ",
      paste(signif(x$ljung_box, 3), collapse = ", "), ")\n", sep = "")
  cat("  correlation of residuals with current and past inputs within 99% band: ",
      x$independent, "\n", sep = "")
  invisible(x)
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_resid <- function(x, ...) {
  ac <- x$acf[x$acf$output == 1, c("lag", "value", "band")]
  cc <- x$ccf[x$ccf$output == 1 & x$ccf$input == 1, c("lag", "value", "band")]
  df <- rbind(data.frame(ac, panel = "residual autocorrelation (output 1)"),
              data.frame(cc, panel = "cor(e[t + lag], u[t]) (output 1, input 1)"))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$lag, y = .data$value)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = .data$band), linetype = "dashed",
                        color = "red") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = -.data$band), linetype = "dashed",
                        color = "red") +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x") +
    ggplot2::labs(x = "lag", y = NULL, title = "Residual tests (99% bands)") +
    ggplot2::theme_bw()
}

#' ARX model order selection
#'
#' Fit ARX models over a grid of orders and report information criteria on
#' the estimation data and, when given, simulation fit on validation data.
#' Prefer validation fit and residual tests to a single criterion; AIC
#' tends to overfit, BIC is more conservative.
#'
#' @param data Estimation data.
#' @param na,nb Candidate orders.
#' @param nk Delay.
#' @param validation Optional validation data.
#' @return A data frame sorted by BIC.
#' @examples
#' set.seed(1)
#' u <- rnorm(400)
#' y <- as.numeric(stats::filter(0.5 * c(0, u[-400]), 0.8, "recursive")) +
#'   rnorm(400, sd = 0.05)
#' parts <- ctrl_idsplit(ctrl_iddata(y, u))
#' head(ctrl_arx_order(parts$estimation, 1:3, 1:3, validation = parts$validation))
#' @export
ctrl_arx_order <- function(data, na = 1:4, nb = 1:4, nk = 1, validation = NULL) {
  grid <- expand.grid(na = na, nb = nb, nk = nk)
  rows <- lapply(seq_len(nrow(grid)), function(r) {
    g <- grid[r, ]
    m <- tryCatch(ctrl_arx(data, g$na, g$nb, g$nk), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    vf <- if (is.null(validation)) NA_real_ else
      ctrl_compare(validation, m)$simulation_fit[1]
    data.frame(na = g$na, nb = g$nb, nk = g$nk, npar = m$npar, aic = m$aic,
               bic = m$bic, validation_fit = vf)
  })
  out <- do.call(rbind, rows)
  out[order(out$bic), , drop = FALSE]
}

#' @export
print.ctrl_idmodel <- function(x, ...) {
  cat("Identified ", toupper(x$method), " model (dt = ", x$dt, ")\n", sep = "")
  if (inherits(x, "ctrl_arx")) {
    cat("  orders na = ", x$na, ", nb = ", x$nb, ", nk = ", x$nk,
        "; noise variance ", signif(x$noise_var, 4), "\n", sep = "")
    cat("  AIC ", signif(x$aic, 6), ", BIC ", signif(x$bic, 6), "\n", sep = "")
  } else {
    cat("  order ", x$order, " (horizon ", x$horizon, ")\n", sep = "")
  }
  print(x$sys)
  invisible(x)
}
