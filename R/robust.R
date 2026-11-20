#' Generalized plants and linear fractional transformations
#'
#' A generalized plant has inputs `(w, u)` (exogenous inputs and controls)
#' and outputs `(z, y)` (performance outputs and measurements):
#' \deqn{\dot x = A x + B_1 w + B_2 u, \quad z = C_1 x + D_{11} w + D_{12} u,
#' \quad y = C_2 x + D_{21} w + D_{22} u.}
#' The last `ncon` inputs are controls and the last `nmeas` outputs are
#' measurements. `ctrl_lft()` closes the loop `u = K y` and returns the map
#' from `w` to `z`; it requires \eqn{I - D_{22} D_K} to be invertible
#' (well-posedness).
#'
#' @param P Generalized plant (`ctrl_ss`).
#' @param K Controller from `y` to `u`.
#' @param nmeas Number of measurements `y`.
#' @param ncon Number of controls `u`.
#' @return `ctrl_lft()` returns the closed loop as a `ctrl_ss` model.
#'   `ctrl_plant_partition()` returns the blocks as a list.
#' @examples
#' P <- ctrl_ss(-1, cbind(1, 1), rbind(1, 1), matrix(0, 2, 2))
#' ctrl_lft(P, ctrl_gain(-2), nmeas = 1, ncon = 1)
#' @export
ctrl_lft <- function(P, K, nmeas, ncon) {
  pp <- ctrl_plant_partition(P, nmeas, ncon)
  k <- as_ss(K)
  # a static gain adopts the plant's time domain, as in ctrl_feedback()
  if (!nrow(k$A)) k$dt <- pp$dt
  if (k$dt != pp$dt) ctrl_abort("`P` and `K` must share the time domain.")
  if (ncol(k$B) != nmeas || nrow(k$C) != ncon) {
    ctrl_abort("`K` must map ", nmeas, " measurements to ", ncon, " controls.")
  }
  with(pp, {
    E1 <- diag(ncon) - k$D %*% D22
    E2 <- diag(nmeas) - D22 %*% k$D
    if (rcond(E1) < 1e-12) ctrl_abort("The interconnection is not well posed (I - DK D22 singular).")
    M <- solve(E1)
    N <- solve(E2)
    Acl <- rbind(cbind(A + B2 %*% M %*% k$D %*% C2, B2 %*% M %*% k$C),
                 cbind(k$B %*% N %*% C2, k$A + k$B %*% N %*% D22 %*% k$C))
    Bcl <- rbind(B1 + B2 %*% M %*% k$D %*% D21, k$B %*% N %*% D21)
    Ccl <- cbind(C1 + D12 %*% M %*% k$D %*% C2, D12 %*% M %*% k$C)
    Dcl <- D11 + D12 %*% M %*% k$D %*% D21
    ctrl_ss(Acl, Bcl, Ccl, Dcl, dt = dt)
  })
}

#' @rdname ctrl_lft
#' @export
ctrl_plant_partition <- function(P, nmeas, ncon) {
  s <- as_ss(P)
  p <- nrow(s$C)
  m <- ncol(s$B)
  nmeas <- as_count(nmeas, "nmeas")
  ncon <- as_count(ncon, "ncon")
  if (nmeas >= p || ncon >= m) {
    ctrl_abort("The plant needs at least one performance output and one ",
               "exogenous input besides the ", nmeas, " measurements and ",
               ncon, " controls.")
  }
  iw <- seq_len(m - ncon)
  iu <- m - ncon + seq_len(ncon)
  iz <- seq_len(p - nmeas)
  iy <- p - nmeas + seq_len(nmeas)
  list(A = s$A, B1 = s$B[, iw, drop = FALSE], B2 = s$B[, iu, drop = FALSE],
       C1 = s$C[iz, , drop = FALSE], C2 = s$C[iy, , drop = FALSE],
       D11 = s$D[iz, iw, drop = FALSE], D12 = s$D[iz, iu, drop = FALSE],
       D21 = s$D[iy, iw, drop = FALSE], D22 = s$D[iy, iu, drop = FALSE],
       dt = s$dt, n = nrow(s$A))
}

synthesis_checks <- function(pp, what) {
  if (pp$dt > 0) {
    ctrl_abort(what, " synthesis is implemented for continuous time plants only.")
  }
  # dimensionless rank of each block: singular values relative to the
  # block's own largest one (independent of the units of z, u, w and y)
  rel_rank <- function(M) {
    d <- svd(M, 0, 0)$d
    if (!length(d) || d[1] == 0) return(0L)
    sum(d > 1e-10 * d[1])
  }
  if (rel_rank(pp$D12) < ncol(pp$D12)) {
    ctrl_abort("D12 must have full column rank (every control must be ",
               "penalized directly in z); this singular problem is not supported.")
  }
  if (rel_rank(pp$D21) < nrow(pp$D21)) {
    ctrl_abort("D21 must have full row rank (every measurement must be ",
               "corrupted by noise in w); this singular problem is not supported.")
  }
  if (!ctrl_controllability(pp$A, pp$B2, discrete = FALSE)$stabilizable) {
    ctrl_abort("(A, B2) is not stabilizable.")
  }
  if (!ctrl_observability(pp$A, pp$C2, discrete = FALSE)$detectable) {
    ctrl_abort("(C2, A) is not detectable.")
  }
}

# Wrap a controller designed for D22 = 0 so it works on the true plant.
d22_shift <- function(K, D22) {
  if (all(D22 == 0)) return(K)
  k <- as_ss(K)
  E <- diag(nrow(k$C)) + k$D %*% D22
  if (rcond(E) < 1e-12) ctrl_abort("Loop shift for D22 is not well posed.")
  Ei <- solve(E)
  # u = K0 (y - D22 u)
  A <- k$A - k$B %*% D22 %*% Ei %*% k$C
  B <- k$B %*% (diag(ncol(k$B)) - D22 %*% Ei %*% k$D)
  C <- Ei %*% k$C
  D <- Ei %*% k$D
  ctrl_ss(A, B, C, D, dt = k$dt)
}

#' H2 optimal output-feedback synthesis
#'
#' Compute the controller minimizing the H2 norm of the closed loop
#' \eqn{T_{zw}} for a continuous time generalized plant (see [ctrl_lft()]).
#'
#' **Assumptions** (checked; violations raise errors): \eqn{D_{11} = 0} (so
#' the closed-loop H2 norm is finite), \eqn{D_{12}} full column rank,
#' \eqn{D_{21}} full row rank, \eqn{(A, B_2)} stabilizable,
#' \eqn{(C_2, A)} detectable, and no invariant zeros of the
#' \eqn{(A, B_2, C_1, D_{12})} and \eqn{(A, B_1, C_2, D_{21})} subsystems on
#' the imaginary axis (detected as Riccati failures).
#'
#' **Algorithm.** Two Riccati equations (Zhou, Doyle and Glover, 1996,
#' Theorem 14.7):
#' \eqn{X_2} = `ctrl_care(A, B2, C1'C1, D12'D12, C1'D12)`,
#' \eqn{F_2 = -(D_{12}^\top D_{12})^{-1}(B_2^\top X_2 + D_{12}^\top C_1)};
#' \eqn{Y_2} = `ctrl_care(A', C2', B1 B1', D21 D21', B1 D21')`,
#' \eqn{L_2 = -(Y_2 C_2^\top + B_1 D_{21}^\top)(D_{21} D_{21}^\top)^{-1}};
#' the controller is \eqn{A_K = A + B_2 F_2 + L_2 C_2}, \eqn{B_K = -L_2},
#' \eqn{C_K = F_2}, \eqn{D_K = 0} (with a loop shift when
#' \eqn{D_{22} \ne 0}). The closed-loop H2 norm is computed independently
#' with [ctrl_norm()] and returned.
#'
#' @param P Generalized plant.
#' @param nmeas,ncon Numbers of measurements and controls.
#' @return An object of class `c("ctrl_synthesis", "ctrl_controller")` with
#'   `K` (controller), `closed_loop`, `norm` (closed-loop H2 norm), `X`,
#'   `Y` and `method`.
#' @references Zhou, K., Doyle, J. C. and Glover, K. (1996). *Robust and
#'   Optimal Control*. Prentice Hall.
#' @examples
#' # LQG as an H2 problem: z = [x1; u], w = [process noise; sensor noise]
#' A <- matrix(c(0, 1, 0, 0), 2, byrow = TRUE)
#' P <- ctrl_ss(A, cbind(diag(2), 0, c(0, 1)),
#'              rbind(c(1, 0), c(0, 0), c(1, 0)),
#'              rbind(c(0, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0)))
#' h2 <- ctrl_h2syn(P, nmeas = 1, ncon = 1)
#' h2$norm
#' @export
ctrl_h2syn <- function(P, nmeas, ncon) {
  pp <- ctrl_plant_partition(P, nmeas, ncon)
  synthesis_checks(pp, "H2")
  if (any(pp$D11 != 0)) {
    ctrl_abort("H2 synthesis needs D11 = 0; otherwise the closed-loop H2 ",
               "norm is infinite.")
  }
  K0 <- with(pp, {
    X <- ctrl_care(A, B2, t(C1) %*% C1, t(D12) %*% D12, t(C1) %*% D12,
                   check_weights = FALSE)$P
    Y <- ctrl_care(t(A), t(C2), B1 %*% t(B1), D21 %*% t(D21), B1 %*% t(D21),
                   check_weights = FALSE)$P
    F2 <- -solve(t(D12) %*% D12, t(B2) %*% X + t(D12) %*% C1)
    L2 <- -t(solve(D21 %*% t(D21), C2 %*% Y + D21 %*% t(B1)))
    structure(ctrl_ss(A + B2 %*% F2 + L2 %*% C2, -L2, F2,
                      matrix(0, ncon, nmeas)), X = X, Y = Y)
  })
  K <- d22_shift(K0, pp$D22)
  cl <- ctrl_lft(P, K, nmeas, ncon)
  if (!ctrl_is_stable(cl)) {
    ctrl_abort("Internal check failed: the H2 controller does not stabilize ",
               "the closed loop.")
  }
  structure(list(K = K, closed_loop = cl, norm = ctrl_norm(cl, "2"),
                 X = attr(K0, "X"), Y = attr(K0, "Y"), method = "h2",
                 nmeas = nmeas, ncon = ncon),
            class = c("ctrl_synthesis", "ctrl_controller"))
}

#' H-infinity output-feedback synthesis (experimental)
#'
#' Compute a stabilizing controller with closed-loop
#' \eqn{\|T_{zw}\|_\infty < \gamma} for a continuous time generalized plant
#' (see [ctrl_lft()]), and optionally search for the smallest achievable
#' \eqn{\gamma} by bisection.
#'
#' **Scope.** Continuous time only, with \eqn{D_{11} = 0},
#' \eqn{D_{12}} full column rank, \eqn{D_{21}} full row rank,
#' \eqn{(A, B_2)} stabilizable and \eqn{(C_2, A)} detectable; arbitrary
#' \eqn{D_{22}} is handled with a loop shift. Plants with
#' \eqn{D_{11} \ne 0} (for example mixed-sensitivity designs with a
#' biproper sensitivity weight) are rejected, as are discrete time plants.
#'
#' **Algorithm.** The plant is first normalized with orthogonal
#' transformations of `z` and `w` and invertible scalings of `u` and `y`
#' so that \eqn{D_{12} = [0; I]} and \eqn{D_{21} = [0, I]} (these leave the
#' H-infinity norm unchanged). For a given \eqn{\gamma}, the two
#' H-infinity Riccati equations of Glover and Doyle (1988) are solved with
#' [ctrl_care()] on indefinite weights; \eqn{\gamma} is feasible when both
#' stabilizing solutions exist, \eqn{X_\infty \succeq 0},
#' \eqn{Y_\infty \succeq 0} and \eqn{\rho(X_\infty Y_\infty) < \gamma^2}. The
#' central controller is then
#' \eqn{A_K = A + \gamma^{-2} B_1 B_1^\top X_\infty + B_2 F_\infty + Z_\infty
#' L_\infty (C_2 + \gamma^{-2} D_{21} B_1^\top X_\infty)},
#' \eqn{B_K = -Z_\infty L_\infty}, \eqn{C_K = F_\infty}, \eqn{D_K = 0}, with
#' \eqn{F_\infty = -(D_{12}^\top C_1 + B_2^\top X_\infty)},
#' \eqn{L_\infty = -(B_1 D_{21}^\top + Y_\infty C_2^\top)} and
#' \eqn{Z_\infty = (I - \gamma^{-2} Y_\infty X_\infty)^{-1}}.
#'
#' **Certification.** Every returned controller is checked independently:
#' the closed loop must be internally stable and its H-infinity norm,
#' computed with [ctrl_norm()] (Bruinsma-Steinbuch, a different
#' algorithm), must be below \eqn{\gamma}. If the check fails, an error is
#' raised instead of returning an uncertified controller. The bisection
#' stops at relative tolerance `tol` and returns the controller for the
#' smallest feasible \eqn{\gamma} found, which is slightly suboptimal by
#' design: the central controller becomes ill-conditioned at the optimum.
#' When the H2 optimal controller already gives a zero closed-loop norm,
#' the optimum zero is attained and that controller is returned with
#' `gamma = 0`.
#'
#' @inheritParams ctrl_h2syn
#' @param gamma Target level; `NULL` searches for the optimum.
#' @param tol Relative bisection tolerance.
#' @param gamma_range Optional `c(lower, upper)` search interval.
#' @return An object of class `c("ctrl_synthesis", "ctrl_controller")` with
#'   `K`, `closed_loop`, `gamma` (the level used), `norm` (certified
#'   closed-loop H-infinity norm), `X`, `Y`, `history` (bisection log) and
#'   `method`.
#' @references Glover, K. and Doyle, J. C. (1988). State-space formulae for
#'   all stabilizing controllers that satisfy an H-infinity norm bound and
#'   relations to risk sensitivity. *Systems and Control Letters*, 11,
#'   167-172.
#' @examples
#' A <- matrix(c(0, 1, 0, 0), 2, byrow = TRUE)
#' P <- ctrl_ss(A, cbind(diag(2), 0, c(0, 1)),
#'              rbind(c(1, 0), c(0, 0), c(1, 0)),
#'              rbind(c(0, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0)))
#' hinf <- ctrl_hinfsyn(P, nmeas = 1, ncon = 1)
#' c(gamma = hinf$gamma, certified_norm = hinf$norm)
#' @export
ctrl_hinfsyn <- function(P, nmeas, ncon, gamma = NULL, tol = 1e-3,
                         gamma_range = NULL) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      tol >= 1) {
    ctrl_abort("`tol` must be a number in (0, 1).")
  }
  pp <- ctrl_plant_partition(P, nmeas, ncon)
  synthesis_checks(pp, "H-infinity")
  if (max(abs(pp$D11)) > 1e-12 * max(1, norm(pp$D12, "2"))) {
    ctrl_abort("This H-infinity implementation requires D11 = 0. General ",
               "D11 (for example a biproper sensitivity weight) is not supported.")
  }
  np <- hinf_normalize(pp)
  feas <- function(g) tryCatch(hinf_central(np, g), error = function(e) NULL)
  history <- data.frame(gamma = numeric(), feasible = logical())
  if (!is.null(gamma)) {
    check_scalar_number(gamma, "gamma", positive = TRUE)
    sol <- feas(gamma)
    history <- rbind(history, data.frame(gamma = gamma, feasible = !is.null(sol)))
    if (is.null(sol)) {
      ctrl_abort("gamma = ", signif(gamma, 6), " is not achievable (a Riccati ",
                 "equation has no stabilizing solution, a solution is not ",
                 "positive semidefinite, or rho(XY) >= gamma^2).")
    }
  } else {
    if (is.null(gamma_range)) {
      h2 <- ctrl_h2syn(P, nmeas, ncon)
      h2inf <- as.numeric(ctrl_norm(h2$closed_loop, "inf"))
      if (h2inf == 0) {
        # the H2 controller already attains the optimum zero: no positive
        # gamma bracket exists, so return it (certified by the norm above)
        return(structure(list(K = h2$K, closed_loop = h2$closed_loop,
                              gamma = 0, norm = 0, X = h2$X, Y = h2$Y,
                              history = history, method = "hinf",
                              nmeas = nmeas, ncon = ncon),
                         class = c("ctrl_synthesis", "ctrl_controller")))
      }
      hi <- 1.01 * h2inf
      lo <- 0
    } else {
      if (!is.numeric(gamma_range) || length(gamma_range) != 2L ||
          any(!is.finite(gamma_range)) || gamma_range[1] < 0 ||
          gamma_range[1] >= gamma_range[2]) {
        ctrl_abort("`gamma_range` must be two finite numbers with ",
                   "0 <= lower < upper.")
      }
      lo <- gamma_range[1]
      hi <- gamma_range[2]
    }
    sol <- feas(hi)
    history <- rbind(history, data.frame(gamma = hi, feasible = !is.null(sol)))
    tries <- 0L
    while (is.null(sol) && tries < 30L) {
      lo <- hi
      hi <- 2 * hi
      sol <- feas(hi)
      history <- rbind(history, data.frame(gamma = hi, feasible = !is.null(sol)))
      tries <- tries + 1L
    }
    if (is.null(sol)) ctrl_abort("No feasible gamma found below ", signif(hi, 4), ".")
    best <- sol
    gbest <- hi
    it <- 0L
    while ((hi - lo) / hi > tol && it < 200L) {
      it <- it + 1L
      mid <- (lo + hi) / 2
      s <- feas(mid)
      history <- rbind(history, data.frame(gamma = mid, feasible = !is.null(s)))
      if (is.null(s)) lo <- mid else {
        hi <- mid
        best <- s
        gbest <- mid
      }
    }
    # back off from the optimum: the central controller is ill-conditioned
    # there and the certification needs a strict margin
    gamma <- gbest * (1 + tol)
    sol <- feas(gamma) %||% best
    if (identical(sol, best)) gamma <- gbest
  }
  K <- d22_shift(hinf_denormalize(sol$K, np), pp$D22)
  cl <- ctrl_lft(P, K, nmeas, ncon)
  if (!ctrl_is_stable(cl)) {
    ctrl_abort("Certification failed: the H-infinity controller does not ",
               "stabilize the closed loop at gamma = ", signif(gamma, 6), ".")
  }
  nrm <- as.numeric(ctrl_norm(cl, "inf"))
  if (!(nrm < gamma * (1 + 1e-8))) {
    ctrl_abort("Certification failed: closed-loop H-infinity norm ",
               signif(nrm, 6), " is not below gamma = ", signif(gamma, 6), ".")
  }
  structure(list(K = K, closed_loop = cl, gamma = gamma, norm = nrm,
                 X = sol$X, Y = sol$Y, history = history, method = "hinf",
                 nmeas = nmeas, ncon = ncon),
            class = c("ctrl_synthesis", "ctrl_controller"))
}

hinf_normalize <- function(pp) {
  p1 <- nrow(pp$D12)
  m2 <- ncol(pp$D12)
  m1 <- ncol(pp$D21)
  p2 <- nrow(pp$D21)
  # z rotation: Theta D12 = [0; Ru]
  q1 <- qr(pp$D12)
  Q1 <- qr.Q(q1, complete = TRUE)
  Ru <- qr.R(q1)[seq_len(m2), seq_len(m2), drop = FALSE]
  Ru <- Ru[, order(q1$pivot), drop = FALSE]
  Theta <- t(Q1[, c(seq_len(p1 - m2) + m2, seq_len(m2)), drop = FALSE])
  # w rotation: D21 Phi = [0, Ry]
  q2 <- qr(t(pp$D21))
  Q2 <- qr.Q(q2, complete = TRUE)
  Ry <- t(qr.R(q2)[seq_len(p2), seq_len(p2), drop = FALSE][, order(q2$pivot),
                                                         drop = FALSE])
  Phi <- Q2[, c(seq_len(m1 - p2) + p2, seq_len(p2)), drop = FALSE]
  Rui <- solve(Ru)
  Ryi <- solve(Ry)
  list(A = pp$A, B1 = pp$B1 %*% Phi, B2 = pp$B2 %*% Rui,
       C1 = Theta %*% pp$C1, C2 = Ryi %*% pp$C2,
       D12 = Theta %*% pp$D12 %*% Rui, D21 = Ryi %*% pp$D21 %*% Phi,
       Rui = Rui, Ryi = Ryi, n = pp$n)
}

hinf_denormalize <- function(K, np) {
  k <- as_ss(K)
  ctrl_ss(k$A, k$B %*% np$Ryi, np$Rui %*% k$C, np$Rui %*% k$D %*% np$Ryi)
}

hinf_central <- function(np, g) {
  A <- np$A
  B1 <- np$B1
  B2 <- np$B2
  C1 <- np$C1
  C2 <- np$C2
  D12 <- np$D12
  D21 <- np$D21
  n <- np$n
  m1 <- ncol(B1)
  m2 <- ncol(B2)
  p1 <- nrow(C1)
  p2 <- nrow(C2)
  B <- cbind(B1, B2)
  R <- blkdiag(-g^2 * diag(m1), diag(m2))
  S <- cbind(matrix(0, n, m1), t(C1) %*% D12)
  X <- ctrl_care(A, B, t(C1) %*% C1, R, S, check_weights = FALSE)$P
  Ct <- rbind(C1, C2)
  Rt <- blkdiag(-g^2 * diag(p1), diag(p2))
  St <- cbind(matrix(0, n, p1), B1 %*% t(D21))
  Y <- ctrl_care(t(A), t(Ct), B1 %*% t(B1), Rt, St, check_weights = FALSE)$P
  tolp <- -1e-9 * max(1, norm(X, "2"), norm(Y, "2"))
  if (min_eig_sym(X) < tolp || min_eig_sym(Y) < tolp) {
    ctrl_abort("Riccati solution not positive semidefinite.")
  }
  rho <- max(Mod(eigen(X %*% Y, only.values = TRUE)$values))
  if (rho >= g^2 * (1 - 1e-9)) ctrl_abort("rho(XY) >= gamma^2.")
  Finf <- -(t(D12) %*% C1 + t(B2) %*% X)
  Linf <- -(B1 %*% t(D21) + Y %*% t(C2))
  Z <- solve(diag(n) - Y %*% X / g^2)
  AK <- A + B1 %*% t(B1) %*% X / g^2 + B2 %*% Finf +
    Z %*% Linf %*% (C2 + D21 %*% t(B1) %*% X / g^2)
  K <- ctrl_ss(AK, -Z %*% Linf, Finf, matrix(0, m2, p2))
  list(K = K, X = X, Y = Y, rho = rho)
}

#' @export
print.ctrl_synthesis <- function(x, ...) {
  if (x$method == "h2") {
    cat("H2 optimal controller: closed-loop H2 norm ", signif(x$norm, 6), "\n",
        sep = "")
  } else {
    cat("H-infinity controller (experimental): gamma = ", signif(x$gamma, 6),
        ", certified closed-loop H-infinity norm ", signif(x$norm, 6), "\n",
        sep = "")
  }
  cat("controller order ", nrow(as_ss(x$K)$A), ", ", x$nmeas,
      " measurement(s), ", x$ncon, " control(s)\n", sep = "")
  invisible(x)
}

#' Mixed-sensitivity generalized plant
#'
#' Build the generalized plant for the mixed-sensitivity problem
#' \deqn{\left\| \begin{pmatrix} W_1 S \\ W_2 K S \\ W_3 T \end{pmatrix}
#' \right\|}
#' with \eqn{S = (I + G K)^{-1}} and \eqn{T = G K S}. The plant has input
#' `(w, u)` with `w` the reference (or output disturbance) and output
#' `(z1, z2, z3, e)` with \eqn{e = w - G u}; the controller uses `e`
#' and the loop convention is \eqn{u = K e}. Any weight may be `NULL`
#' (omitted).
#'
#' For use with [ctrl_hinfsyn()] (which requires \eqn{D_{11} = 0}), `W1`
#' must be strictly proper and `W2` must have full rank feedthrough (so
#' \eqn{D_{12}} has full column rank). `W3 G` must be proper. Plant poles on
#' the imaginary axis (integrators) are not driven by `w` and make the
#' H-infinity Riccati equations singular; shift them slightly into the left
#' half plane for the design and verify the controller on the true plant.
#'
#' @param G Plant.
#' @param W1,W2,W3 Weights on S, KS and T (`ctrl_tf`, `ctrl_ss` or numbers).
#' @return A `ctrl_ss` generalized plant with `nmeas = ncon =` the plant
#'   dimensions.
#' @examples
#' G <- ctrl_tf(1, c(1, 3, 2))
#' W1 <- ctrl_tf(1, c(1, 0.01))
#' P <- ctrl_augw(G, W1 = W1, W2 = 0.1)
#' syn <- ctrl_hinfsyn(P, nmeas = 1, ncon = 1)
#' syn$gamma
#' @export
ctrl_augw <- function(G, W1 = NULL, W2 = NULL, W3 = NULL) {
  g <- as_ss(as_system(G))
  if (g$dt != 0) ctrl_abort("`ctrl_augw()` supports continuous time plants.")
  ny <- nrow(g$C)
  nu <- ncol(g$B)
  ng <- nrow(g$A)
  wt <- function(W, k, name) {
    if (is.null(W)) return(NULL)
    if (is.numeric(W) && !is_ctrl_system(W)) {
      W <- as.matrix(W)
      if (length(W) == 1L) W <- diag(as.numeric(W), k)
      return(ctrl_gain(W))
    }
    w <- as_ss(as_system(W))
    # a static gain has no time domain; a dynamic weight must share the
    # plant's (continuous) one
    if (w$dt != 0) {
      if (nrow(w$A)) {
        ctrl_abort("`", name, "` is a discrete time weight; `ctrl_augw()` needs ",
                   "continuous time weights.")
      }
      w$dt <- 0
    }
    if (ncol(w$B) == 1L && k > 1L) w <- do.call(ctrl_append, rep(list(w), k))
    if (ncol(w$B) != k) ctrl_abort("`", name, "` must have ", k, " inputs.")
    w
  }
  w1 <- wt(W1, ny, "W1")
  w2 <- wt(W2, nu, "W2")
  w3 <- wt(W3, ny, "W3")
  blocks <- Filter(Negate(is.null), list(w1, w2, w3))
  nst <- ng + sum(vapply(blocks, function(b) nrow(b$A), 0L))
  A <- matrix(0, nst, nst)
  Bw <- matrix(0, nst, ny)
  Bu <- matrix(0, nst, nu)
  A[seq_len(ng), seq_len(ng)] <- g$A
  Bu[seq_len(ng), ] <- g$B
  Cz <- NULL
  Dzw <- NULL
  Dzu <- NULL
  off <- ng
  # e = w - Cg xg - Dg u
  Ce <- cbind(-g$C, matrix(0, ny, nst - ng))
  add_block <- function(w, src_C, src_Dw, src_Du) {
    k <- nrow(w$A)
    idx <- off + seq_len(k)
    A[idx, ] <<- A[idx, ] + w$B %*% src_C
    if (k) A[idx, idx] <<- A[idx, idx] + w$A
    Bw[idx, ] <<- w$B %*% src_Dw
    Bu[idx, ] <<- w$B %*% src_Du
    Cblk <- w$D %*% src_C
    if (k) Cblk[, idx] <- Cblk[, idx] + w$C
    Cz <<- rbind(Cz, Cblk)
    Dzw <<- rbind(Dzw, w$D %*% src_Dw)
    Dzu <<- rbind(Dzu, w$D %*% src_Du)
    off <<- off + k
  }
  if (!is.null(w1)) add_block(w1, Ce, diag(ny), -g$D)
  if (!is.null(w2)) add_block(w2, matrix(0, nu, nst), matrix(0, nu, ny), diag(nu))
  if (!is.null(w3)) {
    add_block(w3, cbind(g$C, matrix(0, ny, nst - ng)), matrix(0, ny, ny), g$D)
  }
  if (is.null(Cz)) ctrl_abort("Give at least one weight.")
  C <- rbind(Cz, Ce)
  D <- cbind(rbind(Dzw, diag(ny)), rbind(Dzu, -g$D))
  ctrl_ss(A, cbind(Bw, Bu), C, D)
}

#' Robust stability under unstructured uncertainty
#'
#' Test robust stability of the feedback loop of plant `G` and controller
#' `K` (negative feedback, \eqn{u = -K y}) against a family of
#' unstructured, stable, norm-bounded perturbations \eqn{\Delta} with
#' \eqn{\|\Delta\|_\infty \le 1}:
#'
#' * `"multiplicative_output"`: \eqn{G_\Delta = (I + W \Delta) G}; robustly
#'   stable iff \eqn{\|T W\|_\infty < 1}, \eqn{T = G K (I + G K)^{-1}};
#' * `"multiplicative_input"`: \eqn{G_\Delta = G (I + W \Delta)}; iff
#'   \eqn{\|T_i W\|_\infty < 1}, \eqn{T_i = K G (I + K G)^{-1}};
#' * `"additive"`: \eqn{G_\Delta = G + W \Delta}; iff
#'   \eqn{\|K S W\|_\infty < 1}, \eqn{S = (I + G K)^{-1}}.
#'
#' **Guarantee.** By the small-gain theorem, if the nominal loop is
#' internally stable and the test value is below one, the loop is stable for
#' *every* perturbation in the family (an exact condition for unstructured
#' complex \eqn{\Delta}; the perturbed plants must keep the same number of
#' unstable poles as `G` for the additive and multiplicative forms). The
#' `margin` is the largest scaling of the uncertainty that is tolerated.
#' The result says nothing about structured or real parametric
#' uncertainty, for which it can be conservative, and nothing about robust
#' *performance* (see [ctrl_robust_performance()]).
#'
#' @param G Nominal plant.
#' @param K Controller (negative feedback).
#' @param W Uncertainty weight (system or number).
#' @param type Uncertainty structure.
#' @return A list with `robustly_stable`, `nominally_stable`, `peak` (the
#'   H-infinity norm tested), `margin` (`1 / peak`) and `type`.
#' @examples
#' G <- ctrl_tf(1, c(1, 1))
#' K <- ctrl_tf(10, 1)
#' W <- ctrl_tf(c(1, 0.2), c(0.1, 1))
#' ctrl_robust_stability(G, K, W)
#' @export
ctrl_robust_stability <- function(G, K, W, type = c("multiplicative_output",
                                                    "multiplicative_input",
                                                    "additive")) {
  type <- match.arg(type)
  G <- as_ss(as_system(G))
  K <- as_ss(as_system(K))
  ny <- nrow(G$C)
  nu <- ncol(G$B)
  k <- if (type == "multiplicative_input") nu else ny
  W <- if (is_ctrl_system(W)) as_ss(W) else {
    W <- as.matrix(W)
    if (length(W) == 1L) ctrl_gain(as.numeric(W) * diag(k))
    else if (identical(dim(W), c(k, k))) ctrl_gain(W)
    else ctrl_abort("A numeric `W` must be a scalar or a ", k, " x ", k, " matrix.")
  }
  if (ncol(W$B) == 1L && k > 1L) W <- do.call(ctrl_append, rep(list(W), k))
  L <- ctrl_series(K, G)
  Li <- ctrl_series(G, K)
  cl <- ctrl_feedback(L, ctrl_gain(diag(ny)))
  nominal <- loop_internally_stable(G, K)
  M <- switch(type,
    multiplicative_output = ctrl_series(W, cl),
    multiplicative_input = ctrl_series(W, ctrl_feedback(Li, ctrl_gain(diag(nu)))),
    # K (I + G K)^-1 formed as one feedback interconnection, so the
    # controller's own poles do not remain as hidden modes
    additive = ctrl_series(W, ctrl_feedback(K, G)))
  peak <- if (nominal) as.numeric(ctrl_norm(M, "inf")) else Inf
  list(robustly_stable = nominal && peak < 1, nominally_stable = nominal,
       peak = peak, margin = 1 / peak, type = type)
}

#' Robust performance under multiplicative output uncertainty
#'
#' Test whether the performance specification
#' \eqn{\|W_p S_p\|_\infty < 1} holds for every perturbed plant
#' \eqn{G_p = (I + W_u \Delta) G} with \eqn{\|\Delta\|_\infty \le 1}, where
#' \eqn{S_p = (I + G_p K)^{-1}}. Since
#' \eqn{S_p = S (I + W_u \Delta T)^{-1}}, the condition
#' \deqn{\sup_\omega \left(|W_p(j\omega)|\, \bar\sigma(S(j\omega)) +
#' |W_u(j\omega)|\, \bar\sigma(T(j\omega))\right) < 1}
#' together with nominal internal stability is sufficient. For SISO loops
#' it is also necessary (the classical robust performance condition). For
#' MIMO loops with the scalar weights used here it can be conservative; it
#' is not a structured singular value computation.
#'
#' Each supremum (performance, robustness and their sum) is sampled on a
#' dense frequency grid that includes the endpoints (zero and infinite
#' frequency for continuous models, zero and the Nyquist frequency for
#' discrete models), and the largest local maxima are refined with a
#' one-dimensional search. Sampled peaks are lower bounds on the true
#' suprema, so `robust_performance` is decided on samples. A pass is
#' `certified` when the rigorous upper bound
#' \eqn{\|W_p S\|_\infty + \|W_u T\|_\infty} from [ctrl_norm()] is below
#' one as well; it can be conservative, so an uncertified pass is not a
#' failure.
#'
#' @param G Nominal plant.
#' @param K Controller (negative feedback).
#' @param Wp Scalar performance weight (SISO system or number).
#' @param Wu Scalar uncertainty weight (SISO system or number).
#' @param n Number of grid frequencies.
#' @return A list with `robust_performance` (sampled peak below one),
#'   `certified` (`TRUE` when the rigorous upper bound is also below one),
#'   `upper_bound`, `nominally_stable`, `robust_stability`
#'   (\eqn{\sup |W_u| \bar\sigma(T) < 1}), `nominal_performance`
#'   (\eqn{\sup |W_p| \bar\sigma(S) < 1}), `peak`, `frequency` (where the
#'   peak occurs; `Inf` for the infinite frequency endpoint), `peaks` (the
#'   sampled sum, performance and robustness peaks), `criterion_exact`
#'   (`TRUE` for SISO loops, where the condition is also necessary) and
#'   `curve` (a data frame with `frequency` and `value`).
#' @examples
#' G <- ctrl_tf(1, c(1, 1))
#' K <- ctrl_tf(c(2, 1), c(1, 0))
#' Wp <- ctrl_tf(c(0.5, 0.5), c(1, 0.01))
#' Wu <- ctrl_tf(c(1, 0.2), c(0.1, 1))
#' rp <- ctrl_robust_performance(G, K, Wp, Wu)
#' rp$peak
#' @export
ctrl_robust_performance <- function(G, K, Wp, Wu, n = 2000L) {
  G <- as_ss(as_system(G))
  K <- as_ss(as_system(K))
  ny <- nrow(G$C)
  as_weight <- function(W, name) {
    W <- if (is_ctrl_system(W)) as_system(W) else ctrl_gain(as.numeric(W), dt = G$dt)
    # keep every state of a weight: a near pole-zero cancellation removed by
    # a minimal realization can change the weight at some frequency
    s <- if (is_ctrl_tf(W)) ctrl_tf2ss(W, minimal = FALSE) else as_ss(W)
    if (nrow(s$C) != 1L || ncol(s$B) != 1L) {
      ctrl_abort("`", name, "` must be a scalar (SISO) weight.")
    }
    s
  }
  Wp <- as_weight(Wp, "Wp")
  Wu <- as_weight(Wu, "Wu")
  # one effective clock for the loop and both weights: static parts adopt
  # the time domain of the dynamic ones (a static loop with a discrete
  # weight is evaluated on the unit circle)
  parts <- list(G, K, Wp, Wu)
  dyn <- unique(vapply(Filter(Negate(is_static_model), parts),
                       function(x) x$dt, 0))
  if (length(dyn) > 1L) {
    ctrl_abort("`G`, `K`, `Wp` and `Wu` must share one time domain.")
  }
  clock <- if (length(dyn)) dyn else G$dt
  parts <- lapply(parts, function(x) {
    if (is_static_model(x)) x$dt <- clock
    x
  })
  G <- parts[[1]]
  K <- parts[[2]]
  Wp <- parts[[3]]
  Wu <- parts[[4]]
  nominal <- loop_internally_stable(G, K)
  L <- ctrl_series(K, G)
  S <- ctrl_feedback(ctrl_gain(diag(ny), dt = L$dt), L)
  T <- ctrl_parallel(ctrl_gain(diag(ny), dt = S$dt), -S)
  dt <- S$dt
  smax <- function(R) max(svd(as.matrix(R), 0, 0)$d)
  # values of |Wp| smax(S) and |Wu| smax(T) at frequencies w (Inf allowed
  # for continuous models: the feedthrough terms)
  f <- function(w) {
    fin <- is.finite(w)
    p <- u <- numeric(length(w))
    if (any(fin)) {
      z <- if (dt > 0) exp(1i * w[fin] * dt) else 1i * w[fin]
      RS <- eval_system(S, z)
      RT <- eval_system(T, z)
      p[fin] <- Mod(eval_system(Wp, z)[1, 1, ]) *
        vapply(seq_along(z), function(k) smax(RS[, , k]), 0)
      u[fin] <- Mod(eval_system(Wu, z)[1, 1, ]) *
        vapply(seq_along(z), function(k) smax(RT[, , k]), 0)
    }
    if (any(!fin)) {
      p[!fin] <- abs(Wp$D[1, 1]) * smax(S$D)
      u[!fin] <- abs(Wu$D[1, 1]) * smax(T$D)
    }
    list(p = p, u = u)
  }
  w <- ctrl_freq_grid(ctrl_append(S, Wp, Wu), n = n)
  ends <- if (dt > 0) c(0, pi / dt) else c(0, Inf)
  w <- sort(unique(c(ends, w[w <= (if (dt > 0) pi / dt else Inf)])))
  v <- f(w)
  # refine the largest local maxima of a sampled function
  refine <- function(vals, which) {
    best <- max(vals)
    wbest <- w[which.max(vals)]
    fin <- which(is.finite(w))
    loc <- fin[vapply(fin, function(i) {
      vals[i] >= vals[max(1L, i - 1L)] && vals[i] >= vals[min(length(w), i + 1L)]
    }, logical(1))]
    loc <- loc[order(vals[loc], decreasing = TRUE)][seq_len(min(5L, length(loc)))]
    for (i in loc) {
      lo <- w[max(1L, i - 1L)]
      hi <- w[min(length(w), i + 1L)]
      if (!is.finite(hi)) hi <- 10 * max(w[is.finite(w)])
      if (!(hi > lo)) next
      g <- function(x) {
        r <- f(x)
        switch(which, sum = r$p + r$u, p = r$p, u = r$u)
      }
      o <- stats::optimize(g, c(lo, hi), maximum = TRUE, tol = 1e-10 * max(1, hi))
      if (o$objective > best) {
        best <- o$objective
        wbest <- o$maximum
      }
    }
    c(value = best, frequency = wbest)
  }
  tot <- v$p + v$u
  pk <- refine(tot, "sum")
  pp <- refine(v$p, "p")[["value"]]
  pu <- refine(v$u, "u")[["value"]]
  peak <- pk[["value"]]
  # rigorous upper bound: sup(a + b) <= sup a + sup b, each from the
  # certified H-infinity norm bounds
  upper <- Inf
  if (nominal) {
    nb <- function(W, X) {
      sys <- ctrl_series(X, ctrl_append_n(W, ny))
      if (!ctrl_is_stable(sys)) return(Inf)
      r <- ctrl_norm(sys, "inf")
      attr(r, "bounds")[2] %||% as.numeric(r)
    }
    upper <- tryCatch(nb(Wp, S) + nb(Wu, T), error = function(e) Inf)
    # an upper bound contradicted by a sample is not a bound
    if (upper < peak * (1 - 1e-9)) upper <- Inf
  } else {
    peak <- pp <- pu <- Inf
  }
  list(robust_performance = nominal && peak < 1,
       certified = nominal && upper < 1, upper_bound = upper,
       nominally_stable = nominal,
       robust_stability = nominal && pu < 1,
       nominal_performance = nominal && pp < 1, peak = peak,
       frequency = pk[["frequency"]],
       peaks = c(sum = peak, performance = pp, robustness = pu),
       criterion_exact = ny == 1L && ncol(G$B) == 1L,
       curve = data.frame(frequency = w, value = tot))
}

# Block diagonal of k copies of a SISO model.
ctrl_append_n <- function(W, k) {
  out <- W
  if (k > 1L) for (i in 2:k) out <- ctrl_append(out, W)
  out
}

# Internal stability of the negative feedback loop u = -K y, from the joint
# realization (hidden unstable cancellations are not missed).
loop_internally_stable <- function(G, K) {
  E <- diag(ncol(G$B)) + K$D %*% G$D
  if (rcond(E) < 1e-12) return(FALSE)
  Ei <- solve(E)
  Fx <- -Ei %*% K$D %*% G$C
  Fk <- -Ei %*% K$C
  Acl <- rbind(cbind(G$A + G$B %*% Fx, G$B %*% Fk),
               cbind(K$B %*% (G$C + G$D %*% Fx), K$A + K$B %*% G$D %*% Fk))
  if (!nrow(Acl)) return(TRUE)
  ev <- eigen(Acl, only.values = TRUE)$values
  # a static plant or controller adopts the other's time domain
  if (common_dt(G, K) > 0) all(Mod(ev) < 1) else all(Re(ev) < 0)
}
