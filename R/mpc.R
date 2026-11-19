#' Linear model predictive control
#'
#' Define a constrained linear MPC controller for a discrete time model
#' \eqn{x_{k+1} = A x_k + B u_k + d_k}, \eqn{y_k = C x_k + D u_k}. At each
#' sample the controller solves
#' \deqn{\min \sum_{k=0}^{N-1} \Big[\|x_k - x^r_k\|_Q^2 + \|u_k - u^r_k\|_R^2
#' + \|y_k - y^r_k\|_{Q_y}^2 + \|\Delta u_k\|_{R_d}^2\Big] +
#' \|x_N - x^r_N\|_P^2 + \rho \text{ (slack penalties)}}
#' subject to \eqn{x_0 = x(t)}, the dynamics, input bounds on \eqn{u_k},
#' rate bounds on \eqn{\Delta u_k = u_k - u_{k-1}} (with \eqn{u_{-1}} the
#' previously applied input), state bounds on \eqn{x_1, \ldots, x_N},
#' output bounds on \eqn{y_1, \ldots, y_N} (using \eqn{u_{N-1}} for the
#' feedthrough at \eqn{k = N}) and, when the bounded outputs have direct
#' feedthrough, on \eqn{y_0 = C x_0 + D u_0}, an optional terminal constraint, and
#' \eqn{u_k = u_{N_c - 1}} for \eqn{k \ge N_c} (control horizon). Only
#' \eqn{u_0} is applied (receding horizon).
#'
#' **Indexing.** Decision variables are \eqn{x_0, \ldots, x_N},
#' \eqn{u_0, \ldots, u_{N-1}} and slack variables, in this order. The
#' input \eqn{u_k} acts between \eqn{x_k} and \eqn{x_{k+1}}, so the first
#' predicted state affected by the current decision is \eqn{x_1}.
#'
#' **QP.** The problem is assembled once as a sparse QP in OSQP form
#' \eqn{\min \frac12 z^\top H z + q^\top z}, \eqn{l \le A_c z \le u}; at
#' each sample only `q`, `l` and `u` change, and the OSQP workspace is
#' reused with warm starts. Convexity is verified before the solver is
#' called: `Q`, `P`, `Qy` and `Rd` must be positive semidefinite and `R`
#' positive semidefinite with `R + Rd` positive definite, so the objective
#' is convex and strictly convex in the inputs.
#'
#' **Soft constraints.** With `soft = "states"`, `"outputs"` or
#' `c("states", "outputs")`, the corresponding bounds become
#' \eqn{x_{\min} - s \le x \le x_{\max} + s}, \eqn{s \ge 0}, penalized by
#' \eqn{\rho_q s^2 + \rho_l s}. A large linear weight \eqn{\rho_l} gives an
#' exact penalty: the soft constraints are satisfied whenever the hard
#' problem is feasible. Input and rate bounds are always hard.
#'
#' **Stability.** A successful QP solve does not imply closed-loop
#' stability. Nominal asymptotic stability of the unconstrained-terminal
#' scheme follows when `P` is the LQR cost-to-go (`P = "lqr"`) and the
#' constraints are inactive near the origin; with active constraints one
#' needs a terminal set that is positively invariant under the LQR law,
#' which is not computed automatically. `terminal = "equality"`
#' (\eqn{x_N = x^r_N}) gives recursive feasibility and stability for a
#' constant reference only when \eqn{(x^r, u^r)} is an admissible
#' equilibrium (\eqn{x^r = A x^r + B u^r} with \eqn{u^r} and \eqn{x^r}
#' strictly inside the bounds), the first problem is feasible, and the
#' control horizon equals the horizon with no rate bounds; otherwise the
#' shifted plan need not be feasible.
#'
#' @param plant Discrete time `ctrl_ss` model. Continuous models must be
#'   discretized first with [ctrl_c2d()].
#' @param horizon Prediction horizon `N`.
#' @param Q State weight (scalar, vector or matrix).
#' @param R Input weight (scalar, vector or matrix).
#' @param P Terminal weight: a matrix, `"Q"` (default, same as `Q`) or
#'   `"lqr"`: the DARE solution for the stage cost
#'   \eqn{x^\top Q x + u^\top R u + y^\top Q_y y} with \eqn{y = Cx + Du}
#'   (absolute output weighting, no reference). The move penalty `Rd` is
#'   not part of this cost, so with `Rd` the terminal weight no longer
#'   equals the infinite-horizon cost to go.
#' @param control_horizon Number of free input moves `Nc` (default `N`).
#' @param Rd Weight on input increments.
#' @param Qy Weight on output tracking errors.
#' @param u_min,u_max,du_min,du_max,x_min,x_max,y_min,y_max Bounds (scalars
#'   are recycled; `-Inf`/`Inf` for none).
#' @param soft Which bounds are soft: any of `"states"`, `"outputs"`.
#' @param soft_weight Quadratic and linear slack penalties `c(rho_q, rho_l)`.
#' @param terminal `"none"` or `"equality"`.
#' @param terminal_set Optional list `list(H = , h = )` for the polyhedral
#'   terminal constraint \eqn{H (x_N - x^r_N) \le h}.
#' @param solver_settings List of OSQP settings (see [osqp::osqpSettings()]);
#'   defaults are tight tolerances with polishing.
#' @param accept_inaccurate Accept the OSQP status `"solved inaccurate"`.
#' @return An object of class `c("ctrl_mpc", "ctrl_controller")`.
#' @seealso [ctrl_mpc_step()], [ctrl_simulate_mpc()],
#'   [ctrl_mpc_estimator()], [ctrl_mpc_qp()]
#' @examples
#' plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
#'                           c(1, 0)), 0.1)
#' mpc <- ctrl_mpc(plant, horizon = 20, Q = diag(c(10, 1)), R = 0.1,
#'                 u_min = -1, u_max = 1)
#' step <- ctrl_mpc_step(mpc, x = c(1, 0))
#' step$u
#' sim <- ctrl_simulate_mpc(mpc, x0 = c(1, 0), steps = 60)
#' plot(sim)
#' @export
ctrl_mpc <- function(plant, horizon, Q, R, P = "Q", control_horizon = horizon,
                     Rd = NULL, Qy = NULL, u_min = -Inf, u_max = Inf,
                     du_min = -Inf, du_max = Inf, x_min = -Inf, x_max = Inf,
                     y_min = -Inf, y_max = Inf, soft = character(),
                     soft_weight = c(1e3, 1e4), terminal = c("none", "equality"),
                     terminal_set = NULL, solver_settings = list(),
                     accept_inaccurate = FALSE) {
  terminal <- match.arg(terminal)
  s <- as_ss(plant)
  if (s$dt == 0) {
    ctrl_abort("MPC needs a discrete time model; use ctrl_c2d(plant, dt) first.")
  }
  n <- nrow(s$A)
  m <- ncol(s$B)
  p <- nrow(s$C)
  if (!n) ctrl_abort("MPC needs a model with states.")
  N <- as_count(horizon, "horizon")
  Nc <- as_count(control_horizon, "control_horizon")
  if (Nc > N) ctrl_abort("`control_horizon` cannot exceed `horizon`.")
  Q <- check_psd(expand_weight(Q, n, "Q"), "Q")
  R <- check_psd(expand_weight(R, m, "R"), "R")
  Rd <- if (is.null(Rd)) matrix(0, m, m) else check_psd(expand_weight(Rd, m, "Rd"), "Rd")
  Qy <- if (is.null(Qy)) matrix(0, p, p) else check_psd(expand_weight(Qy, p, "Qy"), "Qy")
  if (min_eig_sym(R + Rd) <= 1e-12 * norm(R + Rd, "2")) {
    ctrl_abort("`R + Rd` must be positive definite so the QP is strictly ",
               "convex in the inputs.")
  }
  if (is.character(P)) {
    P <- match.arg(P, c("Q", "lqr"))
    P <- if (P == "Q") Q else {
      # stage cost including the output term y = Cx + Du
      Qe <- Q + t(s$C) %*% Qy %*% s$C
      Re <- R + t(s$D) %*% Qy %*% s$D
      Ne <- t(s$C) %*% Qy %*% s$D
      if (min_eig_sym(Re) <= 0) ctrl_abort("`P = \"lqr\"` needs positive definite `R + D' Qy D`.")
      ctrl_dare(s$A, s$B, symmetrize(Qe), symmetrize(Re), Ne)$P
    }
  }
  P <- check_psd(expand_weight(P, n, "P"), "P")
  bounds <- list(
    u_min = expand_bound(u_min, m, "u_min"), u_max = expand_bound(u_max, m, "u_max"),
    du_min = expand_bound(du_min, m, "du_min"), du_max = expand_bound(du_max, m, "du_max"),
    x_min = expand_bound(x_min, n, "x_min"), x_max = expand_bound(x_max, n, "x_max"),
    y_min = expand_bound(y_min, p, "y_min"), y_max = expand_bound(y_max, p, "y_max")
  )
  for (nm in c("u", "du", "x", "y")) {
    if (any(bounds[[paste0(nm, "_min")]] > bounds[[paste0(nm, "_max")]])) {
      ctrl_abort("`", nm, "_min` must not exceed `", nm, "_max`.")
    }
  }
  if (!is.numeric(soft_weight) || length(soft_weight) != 2L ||
      any(!is.finite(soft_weight)) || any(soft_weight < 0)) {
    ctrl_abort("`soft_weight` must be two finite nonnegative numbers ",
               "c(rho_q, rho_l).")
  }
  if (length(soft) && !all(soft %in% c("states", "outputs"))) {
    ctrl_abort("`soft` may contain \"states\" and \"outputs\".")
  }
  if (!is.null(terminal_set)) {
    H <- as_real_matrix(terminal_set$H, "terminal_set$H", ncol = n)
    h <- as.numeric(terminal_set$h)
    if (ncol(H) != n || length(h) != nrow(H)) {
      ctrl_abort("`terminal_set` must be list(H = <k x n>, h = <length k>).")
    }
    terminal_set <- list(H = H, h = h)
  }
  settings <- utils::modifyList(list(eps_abs = 1e-9, eps_rel = 1e-9,
                                     max_iter = 50000L, polishing = TRUE,
                                     verbose = FALSE, warm_starting = TRUE),
                                solver_settings)
  mpc <- structure(
    list(plant = s, horizon = N, control_horizon = Nc, Q = Q, R = R, P = P,
         Rd = Rd, Qy = Qy, bounds = bounds, soft = soft,
         soft_weight = soft_weight, terminal = terminal,
         terminal_set = terminal_set, settings = settings,
         accept_inaccurate = accept_inaccurate, dt = s$dt,
         Bd = matrix(0, n, 0), Cd = matrix(0, p, 0)),
    class = c("ctrl_mpc", "ctrl_controller")
  )
  mpc$qp <- mpc_build_qp(mpc)
  mpc$workspace <- new.env(parent = emptyenv())
  mpc
}

as_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != round(x)) {
    ctrl_abort("`", name, "` must be a positive integer.")
  }
  as.integer(x)
}

expand_bound <- function(b, k, name) {
  if (is.null(b)) b <- if (grepl("min", name)) -Inf else Inf
  if (!is.numeric(b) || anyNA(b)) ctrl_abort("`", name, "` must be numeric.")
  if (length(b) == 1L) return(rep(as.numeric(b), k))
  if (length(b) != k) ctrl_abort("`", name, "` must have length 1 or ", k, ".")
  as.numeric(b)
}

# Triplet accumulator for sparse matrices.
trip_new <- function() {
  e <- new.env(parent = emptyenv())
  e$i <- integer()
  e$j <- integer()
  e$x <- numeric()
  e
}

# M is a matrix block, or a vector holding the diagonal of a square block
# (appended directly, without forming the dense matrix).
trip_add <- function(tr, rows, cols, M) {
  if (is.null(dim(M))) {
    keep <- M != 0
    tr$i <- c(tr$i, rows[keep])
    tr$j <- c(tr$j, cols[keep])
    tr$x <- c(tr$x, rep_len(M, length(rows))[keep])
    return(invisible())
  }
  nz <- which(M != 0, arr.ind = TRUE)
  if (!nrow(nz)) return(invisible())
  tr$i <- c(tr$i, rows[nz[, 1]])
  tr$j <- c(tr$j, cols[nz[, 2]])
  tr$x <- c(tr$x, M[nz])
  invisible()
}

mpc_build_qp <- function(mpc) {
  s <- mpc$plant
  A <- s$A
  B <- s$B
  C <- s$C
  D <- s$D
  n <- nrow(A)
  m <- ncol(B)
  p <- nrow(C)
  N <- mpc$horizon
  Nc <- mpc$control_horizon
  b <- mpc$bounds
  xi <- function(k) k * n + seq_len(n)
  nx <- (N + 1L) * n
  ui <- function(k) nx + k * m + seq_len(m)
  nu <- N * m
  # soft rows: one slack per bounded (k, component)
  soft_x <- "states" %in% mpc$soft
  soft_y <- "outputs" %in% mpc$soft
  xb <- which(is.finite(b$x_min) | is.finite(b$x_max))
  yb <- which(is.finite(b$y_min) | is.finite(b$y_max))
  # y_0 = C x_0 + D u_0 depends on the decision u_0 only through D; it is
  # constrained when the bounded outputs have direct feedthrough
  y_stages <- if (length(yb) && any(D[yb, , drop = FALSE] != 0)) 0:N else seq_len(N)
  ns <- (if (soft_x) N * length(xb) else 0L) +
    (if (soft_y) length(y_stages) * length(yb) else 0L)
  si0 <- nx + nu
  nz <- nx + nu + ns
  Hm <- trip_new()
  Ac <- trip_new()
  lo <- numeric()
  hi <- numeric()
  row_tag <- character()
  nrow_c <- 0L
  add_rows <- function(blocks, l, u, tag) {
    k <- length(l)
    rows <- nrow_c + seq_len(k)
    for (bl in blocks) trip_add(Ac, rows, bl$cols, bl$M)
    lo <<- c(lo, l)
    hi <<- c(hi, u)
    row_tag <<- c(row_tag, rep(tag, k))
    nrow_c <<- nrow_c + k
    rows
  }
  # objective Hessian (factor 2 for OSQP's 1/2 z'Hz)
  Qy_x <- t(C) %*% mpc$Qy %*% C
  Qy_xu <- t(C) %*% mpc$Qy %*% D
  Qy_u <- t(D) %*% mpc$Qy %*% D
  for (k in 0:(N - 1L)) {
    trip_add(Hm, xi(k), xi(k), 2 * (mpc$Q + Qy_x))
    trip_add(Hm, ui(k), ui(k), 2 * (mpc$R + Qy_u + mpc$Rd))
    if (any(Qy_xu != 0)) {
      trip_add(Hm, xi(k), ui(k), 2 * Qy_xu)
      trip_add(Hm, ui(k), xi(k), 2 * t(Qy_xu))
    }
    if (k >= 1L) {
      trip_add(Hm, ui(k - 1L), ui(k - 1L), 2 * mpc$Rd)
      trip_add(Hm, ui(k), ui(k - 1L), -2 * mpc$Rd)
      trip_add(Hm, ui(k - 1L), ui(k), -2 * mpc$Rd)
    }
  }
  trip_add(Hm, xi(N), xi(N), 2 * mpc$P)
  if (ns) {
    trip_add(Hm, si0 + seq_len(ns), si0 + seq_len(ns),
             rep(2 * mpc$soft_weight[1], ns))
  }
  # x0 = x(t)
  rows_init <- add_rows(list(list(cols = xi(0), M = diag(n))), rep(0, n),
                        rep(0, n), "init")
  # dynamics x_{k+1} - A x_k - B u_k = d_k
  rows_dyn <- integer()
  for (k in 0:(N - 1L)) {
    r <- add_rows(list(list(cols = xi(k + 1L), M = diag(n)),
                       list(cols = xi(k), M = -A),
                       list(cols = ui(k), M = -B)),
                  rep(0, n), rep(0, n), "dyn")
    rows_dyn <- c(rows_dyn, r)
  }
  # control horizon: u_k = u_{Nc-1}
  if (Nc < N) {
    for (k in Nc:(N - 1L)) {
      add_rows(list(list(cols = ui(k), M = diag(m)),
                    list(cols = ui(Nc - 1L), M = -diag(m))),
               rep(0, m), rep(0, m), "hold")
    }
  }
  # input bounds
  if (any(is.finite(c(b$u_min, b$u_max)))) {
    for (k in 0:(N - 1L)) {
      add_rows(list(list(cols = ui(k), M = diag(m))), b$u_min, b$u_max, "u")
    }
  }
  # rate bounds
  rows_rate0 <- integer()
  if (any(is.finite(c(b$du_min, b$du_max)))) {
    rows_rate0 <- add_rows(list(list(cols = ui(0), M = diag(m))), b$du_min,
                           b$du_max, "rate0")
    if (N > 1L) {
      for (k in 1:(N - 1L)) {
        add_rows(list(list(cols = ui(k), M = diag(m)),
                      list(cols = ui(k - 1L), M = -diag(m))),
                 b$du_min, b$du_max, "rate")
      }
    }
  }
  # state bounds on x_1..x_N
  slack_ptr <- 0L
  if (length(xb)) {
    for (k in 1:N) {
      E <- diag(n)[xb, , drop = FALSE]
      if (soft_x) {
        S <- rep(1, length(xb))
        sc <- si0 + slack_ptr + seq_len(length(xb))
        slack_ptr <- slack_ptr + length(xb)
        add_rows(list(list(cols = xi(k), M = E), list(cols = sc, M = -S)),
                 rep(-Inf, length(xb)), b$x_max[xb], "xs")
        add_rows(list(list(cols = xi(k), M = E), list(cols = sc, M = S)),
                 b$x_min[xb], rep(Inf, length(xb)), "xs")
      } else {
        add_rows(list(list(cols = xi(k), M = E)), b$x_min[xb], b$x_max[xb], "x")
      }
    }
  }
  # output bounds on y_1..y_N (and y_0 with feedthrough)
  rows_y <- integer()
  comp_y <- integer()
  if (length(yb)) {
    for (k in y_stages) {
      Cy <- C[yb, , drop = FALSE]
      Dy <- D[yb, , drop = FALSE]
      ku <- min(k, N - 1L)
      blocks <- list(list(cols = xi(k), M = Cy), list(cols = ui(ku), M = Dy))
      if (soft_y) {
        sc <- si0 + slack_ptr + seq_len(length(yb))
        slack_ptr <- slack_ptr + length(yb)
        S <- rep(1, length(yb))
        r1 <- add_rows(c(blocks, list(list(cols = sc, M = -S))),
                       rep(-Inf, length(yb)), b$y_max[yb], "ys")
        r2 <- add_rows(c(blocks, list(list(cols = sc, M = S))),
                       b$y_min[yb], rep(Inf, length(yb)), "ys")
        rows_y <- c(rows_y, r1, r2)
        comp_y <- c(comp_y, yb, yb)
      } else {
        rows_y <- c(rows_y, add_rows(blocks, b$y_min[yb], b$y_max[yb], "y"))
        comp_y <- c(comp_y, yb)
      }
    }
  }
  if (ns) {
    add_rows(list(list(cols = si0 + seq_len(ns), M = rep(1, ns))), rep(0, ns),
             rep(Inf, ns), "slack")
  }
  rows_term <- integer()
  if (mpc$terminal == "equality") {
    rows_term <- add_rows(list(list(cols = xi(N), M = diag(n))), rep(0, n),
                          rep(0, n), "terminal")
  }
  rows_tset <- integer()
  if (!is.null(mpc$terminal_set)) {
    H <- mpc$terminal_set$H
    rows_tset <- add_rows(list(list(cols = xi(N), M = H)), rep(-Inf, nrow(H)),
                          mpc$terminal_set$h, "terminal_set")
  }
  Hs <- Matrix::sparseMatrix(i = Hm$i, j = Hm$j, x = Hm$x, dims = c(nz, nz))
  Acs <- Matrix::sparseMatrix(i = Ac$i, j = Ac$j, x = Ac$x, dims = c(nrow_c, nz))
  # convexity of the assembled Hessian: each stage block must be PSD
  stage <- rbind(cbind(mpc$Q + Qy_x, Qy_xu), cbind(t(Qy_xu), mpc$R + Qy_u))
  if (min_eig_sym(stage) < -1e-10 * norm(stage, "2")) {
    ctrl_abort("The MPC stage cost is not convex.")
  }
  list(H = Hs, A = Acs, l = lo, u = hi, tags = row_tag, n = n, m = m, p = p,
       N = N, nz = nz, ns = ns, xi = xi, ui = ui, si0 = si0,
       rows_init = rows_init, rows_dyn = rows_dyn, rows_rate0 = rows_rate0,
       rows_y = rows_y, comp_y = comp_y, rows_term = rows_term, rows_tset = rows_tset,
       Qy_x = Qy_x, Qy_xu = Qy_xu, Qy_u = Qy_u)
}

# Expand a reference to N + 1 rows (x), N rows (u, y).
ref_rows <- function(ref, k, width, rows, name) {
  if (is.null(ref)) return(matrix(0, rows, width))
  if (is.function(ref)) {
    out <- t(vapply(k + seq_len(rows) - 1L, function(j) as.numeric(ref(j)),
                    numeric(width)))
    if (width == 1L) out <- matrix(out, ncol = 1)
    return(out)
  }
  if (is.null(dim(ref))) {
    if (length(ref) != width) ctrl_abort("`", name, "` must have length ", width, ".")
    return(matrix(ref, rows, width, byrow = TRUE))
  }
  ref <- as.matrix(ref)
  if (ncol(ref) != width) ctrl_abort("`", name, "` must have ", width, " columns.")
  idx <- pmin(k + seq_len(rows) - 1L, nrow(ref))
  ref[idx, , drop = FALSE]
}

#' Inspect the MPC quadratic program
#'
#' Return the QP data for given parameters, in OSQP form
#' \eqn{\min \frac12 z^\top H z + q^\top z}, \eqn{l \le A z \le u}, together
#' with index helpers. This exposes the exact relation between control
#' variables, predicted states and solver matrices for verification.
#'
#' @param mpc A `ctrl_mpc` object.
#' @param x Current state.
#' @param reference Optional list with `x`, `u`, `y` references (vectors,
#'   matrices with one row per stage, or functions of the stage index).
#' @param u_prev Previously applied input (for rate terms).
#' @param disturbance Known additive disturbance \eqn{d_k} in the dynamics
#'   (vector or matrix with `N` rows).
#' @param k Time index passed to reference functions.
#' @param output_offset Optional known additive output offset (length `p`),
#'   for example an estimated output disturbance: output bounds then apply to
#'   \eqn{C x + D u + \text{offset}}.
#' @return A list with `H`, `q`, `A`, `l`, `u`, `constant` (objective
#'   constant) and the index functions `xi(k)` and `ui(k)`.
#' @examples
#' plant <- ctrl_ss(1, 1, 1, dt = 1)
#' qp <- ctrl_mpc_qp(ctrl_mpc(plant, 3, 1, 1), x = 1)
#' dim(qp$A)
#' @export
ctrl_mpc_qp <- function(mpc, x, reference = NULL, u_prev = NULL,
                        disturbance = NULL, k = 1L, output_offset = NULL) {
  qp <- mpc$qp
  n <- qp$n
  m <- qp$m
  p <- qp$p
  N <- qp$N
  x <- as.numeric(x)
  if (length(x) != n || any(!is.finite(x))) {
    ctrl_abort("`x` must be a finite vector of length ", n, ".")
  }
  u_prev <- if (is.null(u_prev)) rep(0, m) else as.numeric(u_prev)
  if (length(u_prev) != m) ctrl_abort("`u_prev` must have length ", m, ".")
  xr <- ref_rows(reference$x, k, n, N + 1L, "reference$x")
  ur <- ref_rows(reference$u, k, m, N, "reference$u")
  yr <- ref_rows(reference$y, k, p, N, "reference$y")
  d <- if (is.null(disturbance)) matrix(0, N, n) else
    ref_rows(disturbance, 1L, n, N, "disturbance")
  q <- numeric(qp$nz)
  const <- 0
  for (kk in 0:(N - 1L)) {
    xk <- xr[kk + 1L, ]
    uk <- ur[kk + 1L, ]
    yk <- yr[kk + 1L, ]
    q[qp$xi(kk)] <- q[qp$xi(kk)] - 2 * drop(mpc$Q %*% xk) -
      2 * drop(t(mpc$plant$C) %*% mpc$Qy %*% yk)
    q[qp$ui(kk)] <- q[qp$ui(kk)] - 2 * drop(mpc$R %*% uk) -
      2 * drop(t(mpc$plant$D) %*% mpc$Qy %*% yk)
    const <- const + sum(xk * (mpc$Q %*% xk)) + sum(uk * (mpc$R %*% uk)) +
      sum(yk * (mpc$Qy %*% yk))
  }
  xN <- xr[N + 1L, ]
  q[qp$xi(N)] <- q[qp$xi(N)] - 2 * drop(mpc$P %*% xN)
  const <- const + sum(xN * (mpc$P %*% xN))
  q[qp$ui(0)] <- q[qp$ui(0)] - 2 * drop(mpc$Rd %*% u_prev)
  const <- const + sum(u_prev * (mpc$Rd %*% u_prev))
  if (qp$ns) q[qp$si0 + seq_len(qp$ns)] <- mpc$soft_weight[2]
  l <- qp$l
  u <- qp$u
  l[qp$rows_init] <- x
  u[qp$rows_init] <- x
  dv <- as.vector(t(d))
  l[qp$rows_dyn] <- dv
  u[qp$rows_dyn] <- dv
  if (length(qp$rows_rate0)) {
    l[qp$rows_rate0] <- mpc$bounds$du_min + u_prev
    u[qp$rows_rate0] <- mpc$bounds$du_max + u_prev
  }
  if (!is.null(output_offset) && length(qp$rows_y)) {
    off <- as.numeric(output_offset)
    if (length(off) != p) ctrl_abort("`output_offset` must have length ", p, ".")
    # the true output is y = C x + D u + offset: shift the bounds on C x + D u
    l[qp$rows_y] <- l[qp$rows_y] - off[qp$comp_y]
    u[qp$rows_y] <- u[qp$rows_y] - off[qp$comp_y]
  }
  if (length(qp$rows_term)) {
    l[qp$rows_term] <- xN
    u[qp$rows_term] <- xN
  }
  if (length(qp$rows_tset)) {
    u[qp$rows_tset] <- mpc$terminal_set$h + drop(mpc$terminal_set$H %*% xN)
  }
  list(H = qp$H, q = q, A = qp$A, l = l, u = u, constant = const, xi = qp$xi,
       ui = qp$ui, nz = qp$nz, ns = qp$ns, references = list(x = xr, u = ur, y = yr),
       disturbance = d)
}

mpc_solver <- function(mpc, qpd) {
  ws <- mpc$workspace
  if (is.null(ws$model)) {
    settings <- do.call(osqp::osqpSettings, mpc$settings)
    Hu <- methods::as(Matrix::triu(qpd$H), "CsparseMatrix")
    ws$model <- osqp::osqp(P = Hu, q = qpd$q,
                           A = methods::as(qpd$A, "CsparseMatrix"),
                           l = qpd$l, u = qpd$u, pars = settings)
  } else {
    ws$model@Update(q = qpd$q, l = qpd$l, u = qpd$u)
  }
  ws$model
}

#' One MPC step
#'
#' Solve the MPC problem at the current state and return the first input
#' together with the full predicted plan and solver diagnostics. Failed
#' solves (any OSQP status other than `"solved"`, or `"solved inaccurate"`
#' when `accept_inaccurate = TRUE`) are reported with `ok = FALSE` and
#' `u = NA`; the caller decides on a fallback (see [ctrl_simulate_mpc()]).
#'
#' @inheritParams ctrl_mpc_qp
#' @param warm_start Optional primal vector to warm start OSQP (for example
#'   the shifted previous solution).
#' @return A list of class `ctrl_mpc_step` with `u` (input to apply),
#'   `ok`, `status`, `x_plan` (`(N+1) x n`), `u_plan` (`N x m`),
#'   `objective` (including constants), `slack`, `iterations`,
#'   `solve_time`, `primal_residual`, `dual_residual`, `max_violation`
#'   (largest bound violation of the plan, slacks excluded) and
#'   `dynamics_residual`.
#' @examples
#' plant <- ctrl_ss(1, 1, 1, dt = 1)
#' mpc <- ctrl_mpc(plant, 5, 1, 1, u_min = -0.5, u_max = 0.5)
#' ctrl_mpc_step(mpc, x = 3)$u
#' @export
ctrl_mpc_step <- function(mpc, x, reference = NULL, u_prev = NULL,
                          disturbance = NULL, k = 1L, warm_start = NULL,
                          output_offset = NULL) {
  if (!inherits(mpc, "ctrl_mpc")) ctrl_abort("`mpc` must be a `ctrl_mpc` object.")
  qpd <- ctrl_mpc_qp(mpc, x, reference, u_prev, disturbance, k, output_offset)
  model <- mpc_solver(mpc, qpd)
  if (!is.null(warm_start) && length(warm_start) == qpd$nz) {
    model@WarmStart(x = warm_start)
  }
  res <- osqp_solve(model)
  status <- res$info$status
  ok <- identical(status, "solved") ||
    (mpc$accept_inaccurate && identical(status, "solved inaccurate"))
  N <- mpc$horizon
  n <- qpd$nz
  z <- res$x
  plan_ok <- length(z) == qpd$nz && all(is.finite(z))
  xs <- us <- NULL
  obj <- NA_real_
  viol <- NA_real_
  dyn <- NA_real_
  slack <- numeric()
  if (plan_ok) {
    xs <- t(matrix(z[unlist(lapply(0:N, qpd$xi))], nrow = mpc$qp$n))
    us <- t(matrix(z[unlist(lapply(0:(N - 1L), qpd$ui))], nrow = mpc$qp$m))
    obj <- 0.5 * sum(z * as.numeric(qpd$H %*% z)) + sum(qpd$q * z) + qpd$constant
    Az <- drop(as.matrix(qpd$A %*% z))
    hard <- !(mpc$qp$tags %in% c("xs", "ys"))
    viol <- max(0, max((qpd$l - Az)[hard]), max((Az - qpd$u)[hard]))
    d <- qpd$disturbance
    s <- mpc$plant
    dyn <- max(vapply(seq_len(N), function(kk) {
      max(abs(xs[kk + 1L, ] - drop(s$A %*% xs[kk, ] + s$B %*% us[kk, ]) - d[kk, ]))
    }, 0))
    if (qpd$ns) slack <- z[mpc$qp$si0 + seq_len(qpd$ns)]
  }
  u_out <- if (ok && plan_ok) us[1, ] else rep(NA_real_, mpc$qp$m)
  names(u_out) <- mpc$plant$input_names
  out <- list(u = u_out,
              ok = ok && plan_ok, status = status, x_plan = xs, u_plan = us,
              objective = obj, slack = slack, z = z,
              iterations = res$info$iter, solve_time = res$info$solve_time,
              primal_residual = res$info$prim_res,
              dual_residual = res$info$dual_res, max_violation = viol,
              dynamics_residual = dyn)
  class(out) <- "ctrl_mpc_step"
  out
}

#' @export
print.ctrl_mpc_step <- function(x, digits = 4, ...) {
  cat("MPC step: status \"", x$status, "\"", if (!x$ok) " (FAILED)", "\n", sep = "")
  cat("  u = ", paste(signif(x$u, digits), collapse = ", "), "\n", sep = "")
  cat("  objective: ", signif(x$objective, 6), ", iterations: ", x$iterations,
      ", max bound violation: ", signif(x$max_violation, 3), "\n", sep = "")
  invisible(x)
}

#' @export
print.ctrl_mpc <- function(x, ...) {
  b <- x$bounds
  cat("Linear MPC (dt = ", x$dt, "): horizon ", x$horizon, ", control horizon ",
      x$control_horizon, "\n", sep = "")
  cat("  ", x$qp$nz, " decision variables, ", nrow(x$qp$A), " constraint rows\n",
      sep = "")
  show <- function(nm) {
    lo <- b[[paste0(nm, "_min")]]
    hi <- b[[paste0(nm, "_max")]]
    if (any(is.finite(c(lo, hi)))) {
      cat("  ", nm, " in [", paste(lo, collapse = ", "), "] .. [",
          paste(hi, collapse = ", "), "]\n", sep = "")
    }
  }
  for (nm in c("u", "du", "x", "y")) show(nm)
  if (length(x$soft)) cat("  soft: ", paste(x$soft, collapse = ", "), "\n", sep = "")
  if (x$terminal != "none") cat("  terminal constraint: ", x$terminal, "\n", sep = "")
  invisible(x)
}

#' Offset-free MPC estimator
#'
#' Build a disturbance-augmented Kalman filter for offset-free tracking.
#' The plant model is augmented with a constant disturbance
#' \eqn{d_{k+1} = d_k} entering either the outputs (\eqn{y = C x + D u + d},
#' `type = "output"`) or the inputs (\eqn{x_{k+1} = A x + B (u + d)},
#' `type = "input"`). With an observable augmented model and as many
#' disturbances as measured outputs, the closed loop has zero steady-state
#' tracking error for constant references and disturbances whenever it
#' converges, even under model mismatch.
#'
#' At each sample the estimate \eqn{(\hat x, \hat d)} is used to compute a
#' steady-state target \eqn{(x_s, u_s)} from
#' \deqn{\begin{pmatrix} A - I & B \\ C & D \end{pmatrix}
#' \begin{pmatrix} x_s \\ u_s \end{pmatrix} =
#' \begin{pmatrix} -B_d \hat d \\ r - C_d \hat d \end{pmatrix}}
#' (least squares when not square) and the MPC tracks it.
#'
#' @param mpc A `ctrl_mpc` object.
#' @param type `"output"` or `"input"` disturbance model.
#' @param Qn Process noise covariance for the plant states.
#' @param Qd Covariance of the disturbance random walk.
#' @param Rn Measurement noise covariance.
#' @return An object of class `ctrl_mpc_estimator`.
#' @examples
#' plant <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.1)
#' mpc <- ctrl_mpc(plant, 10, Q = 0, Qy = 1, R = 0.01, u_min = -5, u_max = 5)
#' est <- ctrl_mpc_estimator(mpc, "output", Qn = 1e-4, Qd = 1e-3, Rn = 1e-4)
#' @export
ctrl_mpc_estimator <- function(mpc, type = c("output", "input"), Qn, Qd, Rn) {
  type <- match.arg(type)
  s <- mpc$plant
  n <- nrow(s$A)
  m <- ncol(s$B)
  p <- nrow(s$C)
  if (type == "output") {
    nd <- p
    Bd <- matrix(0, n, nd)
    Cd <- diag(p)
  } else {
    nd <- m
    Bd <- s$B
    Cd <- s$D
  }
  Aa <- rbind(cbind(s$A, Bd), cbind(matrix(0, nd, n), diag(nd)))
  Ba <- rbind(s$B, matrix(0, nd, m))
  Ca <- cbind(s$C, Cd)
  aug <- ctrl_ss(Aa, Ba, Ca, s$D, dt = s$dt)
  obs <- ctrl_observability(aug)
  if (!obs$full) {
    ctrl_abort("The disturbance-augmented model is not observable; offset-free ",
               "estimation is impossible with this disturbance model.")
  }
  Qa <- blkdiag(expand_weight(Qn, n, "Qn"), expand_weight(Qd, nd, "Qd"))
  kf <- ctrl_kalman(aug, Qa, Rn)
  structure(list(type = type, Bd = Bd, Cd = Cd, augmented = aug, kalman = kf,
                 n = n, nd = nd),
            class = "ctrl_mpc_estimator")
}

mpc_target <- function(s, Bd, Cd, dhat, r) {
  n <- nrow(s$A)
  m <- ncol(s$B)
  M <- rbind(cbind(s$A - diag(n), s$B), cbind(s$C, s$D))
  rhs <- c(-drop(Bd %*% dhat), r - drop(Cd %*% dhat))
  sol <- if (nrow(M) == ncol(M) && rcond(M) > 1e-12) solve(M, rhs) else
    drop(MASS_ginv(M) %*% rhs)
  list(x = sol[seq_len(n)], u = sol[n + seq_len(m)])
}

#' Closed-loop MPC simulation
#'
#' Run the receding-horizon loop: at each sample, measure (or estimate) the
#' state, solve the MPC problem, apply the first input to the plant and
#' advance. The simulation plant may differ from the prediction model (model
#' mismatch) and may include input disturbances and measurement noise.
#'
#' **Failure handling.** If a QP solve fails, `on_failure` decides what
#' happens; the choice is always recorded in the `fallback` column:
#' `"stop"` (default) raises an error with the solver status, `"hold"`
#' applies the previous input (clipped to the input bounds) and `"shifted"`
#' applies the next move of the last successful plan (then holds its last
#' move). No input is ever applied silently.
#'
#' @param mpc A `ctrl_mpc` object.
#' @param x0 Initial plant state.
#' @param steps Number of samples.
#' @param reference List with `x`, `u` and/or `y` references (see
#'   [ctrl_mpc_qp()]); with an estimator, `y` is the output setpoint.
#' @param plant Simulation plant (default: the prediction model).
#' @param disturbance Input disturbance added to the plant input: a vector
#'   with one value per input (constant), a matrix with one row per sample,
#'   or a function of the sample index. For a single input, a vector with
#'   one value per sample is also accepted. The plant sees the effective
#'   input \eqn{u_k + d_k} in both \eqn{x_{k+1} = A x_k + B(u_k + d_k)} and
#'   \eqn{y_k = C x_k + D(u_k + d_k)}. The measurement used to compute
#'   \eqn{u_k} is taken before \eqn{u_k} is applied, so it sees the input
#'   and disturbance held over the previous sample (zero before the first).
#' @param noise Measurement noise on the outputs (same formats).
#' @param estimator Optional `ctrl_mpc_estimator` for output feedback with
#'   offset-free tracking; without it the true state is measured.
#' @param u0 Input applied before the first sample (for rate constraints).
#' @param on_failure `"stop"`, `"hold"` or `"shifted"`.
#' @param warm_start Warm start each solve with the shifted previous plan.
#' @return An object of class `ctrl_mpc_simulation` with `time`, `x`
#'   (true states, `(steps+1) x n`), `y`, `u`, `xhat`, `dhat`, a `log`
#'   data frame (status, iterations, solve time, objective, max violation,
#'   fallback) and summary counts.
#' @examples
#' plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
#'                           c(1, 0)), 0.1)
#' mpc <- ctrl_mpc(plant, horizon = 20, Q = diag(c(10, 1)), R = 0.1,
#'                 u_min = -1, u_max = 1, x_min = c(-Inf, -0.6),
#'                 x_max = c(Inf, 0.6))
#' sim <- ctrl_simulate_mpc(mpc, x0 = c(2, 0), steps = 80)
#' summary(sim)
#' @export
ctrl_simulate_mpc <- function(mpc, x0, steps, reference = NULL, plant = NULL,
                              disturbance = NULL, noise = NULL, estimator = NULL,
                              u0 = NULL, on_failure = c("stop", "hold", "shifted"),
                              warm_start = TRUE) {
  on_failure <- match.arg(on_failure)
  if (!inherits(mpc, "ctrl_mpc")) ctrl_abort("`mpc` must be a `ctrl_mpc` object.")
  model <- mpc$plant
  sim <- if (is.null(plant)) model else as_ss(plant)
  if (sim$dt != model$dt) ctrl_abort("The simulation plant must share the MPC sample time.")
  if (ncol(sim$B) != ncol(model$B) || nrow(sim$C) != nrow(model$C)) {
    ctrl_abort("The simulation plant must have the inputs and outputs of the model.")
  }
  steps <- as_count(steps, "steps")
  n <- nrow(model$A)
  ns <- nrow(sim$A)
  m <- ncol(model$B)
  p <- nrow(model$C)
  N <- mpc$horizon
  x <- as.numeric(x0)
  if (length(x) != ns) ctrl_abort("`x0` must have length ", ns, ".")
  if (is.null(estimator) && ns != n) {
    ctrl_abort("Without an estimator the simulation plant must have the model's states.")
  }
  sig <- function(v, k, width) {
    if (is.null(v)) return(rep(0, width))
    if (is.function(v)) return(as.numeric(v(k)))
    if (is.null(dim(v))) {
      v <- as.numeric(v)
      if (length(v) == width) return(v)
      if (length(v) == 1L) return(rep(v, width))
      # a single channel given as a sequence over the samples
      if (width == 1L && length(v) >= steps) return(v[k])
      ctrl_abort("A signal vector must have length 1 or ", width,
                 " (one value per channel); give a sequence over the samples ",
                 "as a matrix with one row per sample.")
    }
    if (ncol(v) != width) ctrl_abort("A signal matrix must have ", width, " columns.")
    if (nrow(v) < steps) ctrl_abort("A signal matrix must have at least ", steps, " rows.")
    v[k, ]
  }
  X <- matrix(0, steps + 1L, ns)
  Y <- matrix(0, steps, p)
  U <- matrix(0, steps, m)
  Xh <- matrix(NA_real_, steps, n)
  Dh <- NULL
  X[1, ] <- x
  u_prev <- if (is.null(u0)) rep(0, m) else as.numeric(u0)
  log <- data.frame(k = seq_len(steps), status = "", ok = FALSE, iterations = 0L,
                    solve_time = 0, objective = NA_real_, max_violation = NA_real_,
                    max_slack = 0, fallback = "none", stringsAsFactors = FALSE)
  last_plan <- NULL
  plan_age <- 0L
  zprev <- NULL
  # estimator state
  if (!is.null(estimator)) {
    kf <- estimator$kalman
    xa <- rep(0, n + estimator$nd)
    Dh <- matrix(NA_real_, steps, estimator$nd)
  }
  d_prev <- rep(0, m)
  for (k in seq_len(steps)) {
    # the measurement is taken before the new input is applied, so it sees
    # the input and disturbance held over the previous sample
    dk <- sig(disturbance, k, m)
    y <- drop(sim$C %*% x + sim$D %*% (u_prev + d_prev)) + sig(noise, k, p)
    ref <- reference
    dist_model <- NULL
    y_off <- NULL
    if (!is.null(estimator)) {
      # measurement update of the augmented filter (D term uses u_prev as
      # the held input until the new input is computed)
      e <- y - drop(estimator$augmented$C %*% xa + estimator$augmented$D %*% u_prev)
      xa <- xa + drop(kf$M %*% e)
      xhat <- xa[seq_len(n)]
      dhat <- xa[n + seq_len(estimator$nd)]
      Dh[k, ] <- dhat
      r_y <- ref_rows(reference$y, k, p, 1L, "reference$y")[1, ]
      tg <- mpc_target(model, estimator$Bd, estimator$Cd, dhat, r_y)
      ref <- list(x = tg$x, u = tg$u, y = r_y - drop(estimator$Cd %*% dhat))
      dist_model <- matrix(drop(estimator$Bd %*% dhat), N, n, byrow = TRUE)
      y_off <- drop(estimator$Cd %*% dhat)
    } else {
      xhat <- x
    }
    Xh[k, ] <- xhat
    ws <- NULL
    if (warm_start && !is.null(zprev)) ws <- shift_plan(zprev, mpc)
    st <- ctrl_mpc_step(mpc, xhat, ref, u_prev, dist_model, k, ws, y_off)
    log$status[k] <- st$status
    log$ok[k] <- st$ok
    log$iterations[k] <- as.integer(st$iterations %||% NA_integer_)
    log$solve_time[k] <- st$solve_time %||% NA_real_
    log$objective[k] <- st$objective
    log$max_violation[k] <- st$max_violation
    log$max_slack[k] <- if (length(st$slack)) max(st$slack) else 0
    if (st$ok) {
      u <- st$u
      last_plan <- st$u_plan
      plan_age <- 0L
      zprev <- st$z
    } else {
      if (on_failure == "stop") {
        ctrl_abort("MPC solve failed at sample ", k, " with OSQP status \"",
                   st$status, "\". Use `on_failure = \"hold\"` or \"shifted\" ",
                   "to continue with an explicit fallback.")
      }
      plan_age <- plan_age + 1L
      if (on_failure == "shifted" && !is.null(last_plan)) {
        u <- last_plan[min(plan_age + 1L, nrow(last_plan)), ]
        log$fallback[k] <- "shifted"
      } else {
        u <- pmin(pmax(u_prev, mpc$bounds$u_min), mpc$bounds$u_max)
        log$fallback[k] <- "hold"
      }
      zprev <- NULL
    }
    U[k, ] <- u
    Y[k, ] <- drop(sim$C %*% x + sim$D %*% (u + dk))
    x <- drop(sim$A %*% x + sim$B %*% (u + dk))
    X[k + 1L, ] <- x
    if (!is.null(estimator)) {
      xa <- drop(estimator$augmented$A %*% xa + estimator$augmented$B %*% u)
    }
    u_prev <- u
    d_prev <- dk
  }
  colnames(X) <- sim$state_names
  colnames(U) <- model$input_names
  colnames(Y) <- model$output_names
  structure(list(time = (0:steps) * model$dt, x = X, y = Y, u = U, xhat = Xh,
                 dhat = Dh, log = log, bounds = mpc$bounds,
                 failures = sum(!log$ok), dt = model$dt,
                 reference = reference),
            class = "ctrl_mpc_simulation")
}

shift_plan <- function(z, mpc) {
  qp <- mpc$qp
  N <- qp$N
  out <- z
  for (k in 0:(N - 1L)) out[qp$xi(k)] <- z[qp$xi(k + 1L)]
  if (N > 1L) for (k in 0:(N - 2L)) out[qp$ui(k)] <- z[qp$ui(k + 1L)]
  out
}

#' @export
summary.ctrl_mpc_simulation <- function(object, ...) {
  lg <- object$log
  b <- object$bounds
  U <- object$u
  X <- object$x[-1, , drop = FALSE]
  vu <- max(0, max(sweep(U, 2, b$u_max)), max(sweep(-U, 2, -b$u_min)))
  vx <- if (ncol(X) == length(b$x_max)) {
    max(0, max(sweep(X, 2, b$x_max)), max(sweep(-X, 2, -b$x_min)))
  } else NA_real_
  out <- list(steps = nrow(lg), failures = sum(!lg$ok),
              failure_rate = mean(!lg$ok),
              fallbacks = table(lg$fallback),
              statuses = table(lg$status),
              mean_solve_time = mean(lg$solve_time, na.rm = TRUE),
              max_solve_time = max(lg$solve_time, na.rm = TRUE),
              mean_iterations = mean(lg$iterations, na.rm = TRUE),
              max_input_violation = vu, max_state_violation = vx,
              max_slack = max(lg$max_slack))
  class(out) <- "summary.ctrl_mpc_simulation"
  out
}

#' @export
print.summary.ctrl_mpc_simulation <- function(x, ...) {
  cat("MPC simulation: ", x$steps, " samples, ", x$failures, " failed solves (",
      signif(100 * x$failure_rate, 3), "%)\n", sep = "")
  cat("  solver statuses: ", paste(names(x$statuses), x$statuses, sep = " = ",
                                   collapse = ", "), "\n", sep = "")
  cat("  solve time: mean ", signif(x$mean_solve_time, 3), " s, max ",
      signif(x$max_solve_time, 3), " s; mean iterations ",
      signif(x$mean_iterations, 4), "\n", sep = "")
  cat("  max input bound violation: ", signif(x$max_input_violation, 3),
      "\n  max state bound violation (true plant): ",
      signif(x$max_state_violation, 3), "\n  max slack: ",
      signif(x$max_slack, 3), "\n", sep = "")
  invisible(x)
}

#' @export
print.ctrl_mpc_simulation <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

#' @rdname ctrl_plots
#' @export
plot.ctrl_mpc_simulation <- function(x, ...) {
  steps <- nrow(x$u)
  t <- x$time
  rows <- list()
  for (j in seq_len(ncol(x$x))) {
    rows[[length(rows) + 1L]] <- data.frame(time = t, value = x$x[, j],
                                            panel = paste("state", colnames(x$x)[j]),
                                            kind = "trajectory")
  }
  for (j in seq_len(ncol(x$u))) {
    rows[[length(rows) + 1L]] <- data.frame(time = t[seq_len(steps)],
                                            value = x$u[, j],
                                            panel = paste("input", colnames(x$u)[j]),
                                            kind = "trajectory")
  }
  df <- do.call(rbind, rows)
  lim <- list()
  b <- x$bounds
  add_lim <- function(vals, prefix, names) {
    for (j in seq_along(vals)) if (is.finite(vals[j])) {
      lim[[length(lim) + 1L]] <<- data.frame(panel = paste(prefix, names[j]),
                                             value = vals[j])
    }
  }
  add_lim(b$u_min, "input", colnames(x$u))
  add_lim(b$u_max, "input", colnames(x$u))
  if (length(b$x_min) == ncol(x$x)) {
    add_lim(b$x_min, "state", colnames(x$x))
    add_lim(b$x_max, "state", colnames(x$x))
  }
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_step()
  if (length(lim)) {
    g <- g + ggplot2::geom_hline(data = do.call(rbind, lim),
                                 ggplot2::aes(yintercept = .data$value),
                                 linetype = "dashed", color = "red")
  }
  fails <- x$log$k[!x$log$ok]
  if (length(fails)) {
    g <- g + ggplot2::geom_vline(xintercept = t[fails], color = "orange",
                                 alpha = 0.5)
  }
  g + ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = time_label(x$dt), y = NULL,
                  title = "MPC closed loop (dashed: bounds)") +
    ggplot2::theme_bw()
}
