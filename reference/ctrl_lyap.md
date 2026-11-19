# Lyapunov and Stein equations

Solve the continuous time Lyapunov equation \$\$A X + X A^\top + Q =
0\$\$ or the discrete time Lyapunov (Stein) equation \$\$A X A^\top -
X + Q = 0.\$\$

## Usage

``` r
ctrl_lyap(A, Q, discrete = FALSE, tol = 1e-12)
```

## Arguments

- A:

  Square matrix.

- Q:

  Square matrix of the same size (usually symmetric).

- discrete:

  `FALSE` (default) for the continuous equation, `TRUE` for the discrete
  equation.

- tol:

  Relative uniqueness tolerance.

## Value

The solution matrix `X`, with attributes `residual` (relative Frobenius
residual of the equation) and `separation` (the uniqueness measure).
When `Q` is symmetric the solution is symmetrized.

## Details

Both solvers reduce `A` to real Schur form with LAPACK. The continuous
equation is then solved with the Bartels-Stewart algorithm (LAPACK
`dtrsyl`); the discrete equation with a block back substitution over the
quasi-triangular Schur factor (Barraud's method). No Kronecker
vectorization of the full problem is used, so the cost is \\O(n^3)\\.

The solution is unique exactly when \\\lambda_i + \lambda_j \ne 0\\
(continuous) or \\\lambda_i \lambda_j \ne 1\\ (discrete) for all pairs
of eigenvalues of `A`. The separation is \\\min \|\lambda_i +
\lambda_j\| / (\|\lambda_i\| + \|\lambda_j\|)\\ (continuous; zero when a
sum is within rounding of the spectral radius) or \\\min \|1 - \lambda_i
\lambda_j\| / \max(1, \max\|\lambda\|)^2\\ (discrete). If it is below
`tol`, an error is raised instead of returning an unreliable solution.

## Examples

``` r
A <- matrix(c(-1, 2, 0, -3), 2)
X <- ctrl_lyap(A, diag(2))
A %*% X + X %*% t(A) + diag(2)  # close to zero
#>      [,1]         [,2]
#> [1,]    0 0.000000e+00
#> [2,]    0 1.110223e-16
Xd <- ctrl_lyap(A / 4, diag(2), discrete = TRUE)
attr(Xd, "residual")
#> [1] 0
```
