test_that("continuous and discrete Lyapunov solvers match Kronecker references", {
  set.seed(2)
  n <- 7
  A <- matrix(rnorm(n * n), n)
  Ac <- A - (max(Re(eigen(A)$values)) + 1) * diag(n)
  Q <- crossprod(matrix(rnorm(n * n), n))
  X <- ctrl_lyap(Ac, Q)
  ref <- matrix(solve(kronecker(diag(n), Ac) + kronecker(Ac, diag(n)), -as.vector(Q)), n)
  expect_equal(X, ref, tolerance = 1e-9, ignore_attr = TRUE)
  expect_lt(attr(X, "residual"), 1e-12)
  Ad <- A / (1.1 * max(Mod(eigen(A)$values)))
  Xd <- ctrl_lyap(Ad, Q, discrete = TRUE)
  refd <- matrix(solve(diag(n^2) - kronecker(Ad, Ad), as.vector(Q)), n)
  expect_equal(Xd, refd, tolerance = 1e-9, ignore_attr = TRUE)
  expect_lt(attr(Xd, "residual"), 1e-12)
  expect_equal(Xd, t(Xd), ignore_attr = TRUE)
})

test_that("Lyapunov uniqueness failures raise errors", {
  expect_error(ctrl_lyap(matrix(c(0, 1, -1, 0), 2), diag(2)), "unique")
  expect_error(ctrl_lyap(diag(c(2, 0.5)), diag(2), discrete = TRUE), "unique")
  expect_error(ctrl_lyap(diag(2), diag(3)), "2 x 2")
})

test_that("Sylvester and Gramians", {
  A <- matrix(c(-1, 1, 0, -2), 2)
  B <- matrix(-3)
  C <- matrix(c(1, 2), 2)
  X <- ctrl_sylvester(A, B, C)
  expect_equal(A %*% X + X %*% B, C, tolerance = 1e-12, ignore_attr = TRUE)
  sys <- ctrl_ss(-1, 1, 1)
  expect_equal(drop(ctrl_gram(sys)), 0.5, tolerance = 1e-14, ignore_attr = TRUE)
  expect_error(ctrl_gram(ctrl_ss(1, 1, 1)), "stable")
})
