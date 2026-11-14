# Timing benchmarks for controlcompR. Not run on CRAN or in R CMD check.
#   Rscript benchmarks/performance.R
library(controlcompR)
set.seed(1)
time_it <- function(expr, reps = 5L) {
  f <- substitute(expr)
  env <- parent.frame()
  t <- vapply(seq_len(reps), function(i) system.time(eval(f, env))[["elapsed"]], 0)
  c(median = stats::median(t), max = max(t))
}
rows <- list()
for (n in c(10, 50, 100)) {
  A <- matrix(rnorm(n * n), n) / sqrt(n) - 1.5 * diag(n)
  B <- matrix(rnorm(n * 2), n)
  Q <- diag(n)
  R <- diag(2)
  rows[[length(rows) + 1L]] <- c(task = paste("care n =", n),
                                 time_it(ctrl_care(A, B, Q, R)))
  rows[[length(rows) + 1L]] <- c(task = paste("lyap n =", n),
                                 time_it(ctrl_lyap(A, Q)))
}
plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                          c(1, 0)), 0.1)
for (N in c(10, 30, 100)) {
  mpc <- ctrl_mpc(plant, N, Q = diag(c(10, 1)), R = 0.1, u_min = -1, u_max = 1,
                  x_min = c(-Inf, -0.6), x_max = c(Inf, 0.6))
  sim <- ctrl_simulate_mpc(mpc, c(2, 0), 100)
  s <- summary(sim)
  rows[[length(rows) + 1L]] <- c(task = paste("MPC solve, N =", N),
                                 median = stats::median(sim$log$solve_time),
                                 max = s$max_solve_time)
}
out <- as.data.frame(do.call(rbind, rows))
print(out, row.names = FALSE)
cat("Times in seconds; MPC times are OSQP solve times per sample.\n")
