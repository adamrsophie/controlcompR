mass_spring <- function() {
  ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
}

random_stable <- function(n, m, p, seed = 1, dt = 0) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), n)
  if (dt > 0) {
    A <- A / (1.2 * max(Mod(eigen(A)$values)))
  } else {
    A <- A - (max(Re(eigen(A)$values)) + 0.5) * diag(n)
  }
  ctrl_ss(A, matrix(rnorm(n * m), n), matrix(rnorm(p * n), p),
          matrix(0, p, m), dt = dt)
}

sort_complex_test <- function(z) z[order(Re(z), Im(z))]
