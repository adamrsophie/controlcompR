# Timing and memory benchmarks for controlcompR. Not run on CRAN or in
# R CMD check.
#   Rscript benchmarks/performance.R
#
# Each row reports the median and maximum elapsed time over `reps` runs and
# the peak additional R heap memory (Mb) seen by the garbage collector
# during one run. Peak memory comes from gc(reset = TRUE), is approximate,
# and excludes native allocations (LAPACK, OSQP and the package's C code);
# it is not the process resident set size.
library(controlcompR)
set.seed(1)

measure <- function(expr, reps = 5L) {
  f <- substitute(expr)
  env <- parent.frame()
  t <- vapply(seq_len(reps), function(i) system.time(eval(f, env))[["elapsed"]], 0)
  g0 <- gc(reset = TRUE)
  eval(f, env)
  g1 <- gc()
  c(median = stats::median(t), max = max(t),
    peak_mb = max(0, sum(g1[, 6]) - sum(g0[, 2])))
}
rows <- list()
add <- function(task, m) rows[[length(rows) + 1L]] <<- c(task = task, signif(m, 4))

random_stable <- function(n, m = 2L, p = 2L) {
  A <- matrix(rnorm(n * n), n) / sqrt(n) - 1.5 * diag(n)
  list(A = A, B = matrix(rnorm(n * m), n), C = matrix(rnorm(p * n), p))
}

for (n in c(10, 50, 100)) {
  M <- random_stable(n)
  # model construction and validation
  add(paste("ctrl_ss construction, n =", n), measure(ctrl_ss(M$A, M$B, M$C)))
  sys <- ctrl_ss(M$A, M$B, M$C)
  add(paste("care, n =", n), measure(ctrl_care(M$A, M$B, diag(n), diag(2))))
  add(paste("lyap, n =", n), measure(ctrl_lyap(M$A, diag(n))))
  add(paste("H-infinity norm, n =", n), measure(ctrl_norm(sys, "inf")))
  # time simulation
  add(paste("step response, n =", n), measure(ctrl_step(sys)))
  t <- seq(0, 10, by = 0.01)
  add(paste("simulate 1001 samples, n =", n),
      measure(ctrl_simulate(sys, input = cbind(sin(t), cos(t)), times = t)))
  # Kalman filtering
  sd <- ctrl_c2d(sys, 0.05)
  y <- matrix(rnorm(1000 * 2), 1000)
  add(paste("Kalman filter 1000 samples, n =", n),
      measure(ctrl_kalman_filter(sd, y, Qn = diag(n) * 1e-3, Rn = diag(2) * 0.01),
              reps = 3L))
}

plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                          c(1, 0)), 0.1)
for (N in c(10, 30, 100)) {
  # QP assembly (constructor builds the sparse QP once)
  build <- function() {
    ctrl_mpc(plant, N, Q = diag(c(10, 1)), R = 0.1, u_min = -1, u_max = 1,
             x_min = c(-Inf, -0.6), x_max = c(Inf, 0.6))
  }
  add(paste("MPC QP assembly, N =", N), measure(build()))
  mpc <- build()
  # parameter update: QP vectors for a new state and reference
  add(paste("MPC parameter update, N =", N),
      measure(ctrl_mpc_qp(mpc, c(1, 0), reference = list(x = c(0.5, 0))), reps = 20L))
  # first solve (creates the OSQP workspace) and later solves (update only)
  first <- system.time(ctrl_mpc_step(mpc, c(2, 0)))[["elapsed"]]
  add(paste("MPC first step (workspace setup), N =", N),
      c(median = first, max = first, peak_mb = NA))
  add(paste("MPC later step (workspace reuse), N =", N),
      measure(ctrl_mpc_step(mpc, c(1.5, 0.1)), reps = 20L))
  sim <- ctrl_simulate_mpc(mpc, c(2, 0), 100)
  add(paste("MPC OSQP solve time per sample, N =", N),
      c(median = stats::median(sim$log$solve_time),
        max = summary(sim)$max_solve_time, peak_mb = NA))
}

# moving horizon estimation with and without warm starts
sys <- ctrl_ss(matrix(c(1, 0.1, 0, 0.95), 2, byrow = TRUE), c(0, 0.1), c(1, 0), dt = 1)
y <- cumsum(rnorm(200, sd = 0.1))
mhe <- ctrl_mhe(sys, 10, Qn = 1e-4, Rn = 0.01, x_min = c(-5, -5), x_max = c(5, 5),
                solver_settings = list(check_termination = 1L))
for (ws in c(FALSE, TRUE)) {
  add(paste("MHE 200 samples, warm start =", ws),
      measure(ctrl_mhe_run(mhe, y, warm_start = ws), reps = 3L))
}

out <- as.data.frame(do.call(rbind, rows))
print(out, row.names = FALSE)
cat("Times in seconds; peak_mb is the peak extra R heap memory during one run",
    "(from gc(); native allocations by LAPACK, OSQP or C code are not counted).\n")
