test_that("scalar CARE matches sqrt(2) - 1", {
  r <- ctrl_care(-1, 1, 1, 1)
  expect_equal(drop(r$P), sqrt(2) - 1, tolerance = 1e-12)
  expect_equal(drop(r$K), sqrt(2) - 1, tolerance = 1e-12)
  expect_lt(r$residual, 1e-12)
})

test_that("scalar DARE matches the golden ratio", {
  r <- ctrl_dare(1, 1, 1, 1)
  P <- (1 + sqrt(5)) / 2
  expect_equal(drop(r$P), P, tolerance = 1e-12)
  expect_equal(drop(r$K), P / (1 + P), tolerance = 1e-12)
  expect_lt(Mod(r$poles), 1)
})

test_that("CARE with cross term matches the closed form", {
  r <- ctrl_care(-1, 1, 1, 1, S = 0.5)
  expect_equal(drop(r$P), (-3 + sqrt(12)) / 2, tolerance = 1e-12)
})

test_that("random CARE and DARE solutions satisfy the equations", {
  for (seed in 1:5) {
    set.seed(seed)
    n <- 6
    m <- 2
    A <- matrix(rnorm(n * n), n)
    B <- matrix(rnorm(n * m), n)
    Q <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
    R <- diag(m) + 0.1
    S <- matrix(rnorm(n * m), n) * 0.1
    rc <- ctrl_care(A, B, Q, R, S)
    expect_lt(rc$residual, 1e-10)
    expect_true(all(Re(rc$poles) < 0))
    rd <- ctrl_dare(A, B, Q, R, S)
    expect_lt(rd$residual, 1e-10)
    expect_true(all(Mod(rd$poles) < 1))
    # Hamiltonian eigenvector reference for the CARE without cross term
    rc0 <- ctrl_care(A, B, Q, R)
    H <- rbind(cbind(A, -B %*% solve(R, t(B))), cbind(-Q, -t(A)))
    e <- eigen(H)
    V <- e$vectors[, Re(e$values) < 0]
    Pref <- Re(V[(n + 1):(2 * n), ] %*% solve(V[1:n, ]))
    expect_equal(rc0$P, Pref, tolerance = 1e-8)
  }
})

test_that("DARE works with singular A", {
  A <- matrix(c(0, 1, 0, 0), 2)
  r <- ctrl_dare(A, c(1, 0), diag(2), 1)
  expect_lt(r$residual, 1e-12)
})

test_that("Riccati prerequisites are checked", {
  expect_error(ctrl_care(-1, 1, 1, 0), "positive definite")
  expect_error(ctrl_care(-1, 1, -1, 1), "positive semidefinite")
  # unstabilizable: unstable mode not reachable
  expect_error(ctrl_care(diag(c(1, -1)), c(0, 1), diag(2), 1), "stabilizing")
  # unobservable mode on the imaginary axis
  expect_error(ctrl_care(matrix(c(0, 1, -1, 0), 2), c(0, 1), diag(c(0, 0)), 1),
               "stabilizing")
})
