# Algebraic Riccati equations

`ctrl_care()` computes the stabilizing solution of the continuous time
algebraic Riccati equation \$\$A^\top P + P A - (P B + S) R^{-1} (B^\top
P + S^\top) + Q = 0,\$\$ and `ctrl_dare()` the stabilizing solution of
the discrete time equation \$\$A^\top P A - P - (A^\top P B + S)(R +
B^\top P B)^{-1} (B^\top P A + S^\top) + Q = 0.\$\$

## Usage

``` r
ctrl_care(A, B, Q, R, S = NULL, check_weights = TRUE)

ctrl_dare(A, B, Q, R, S = NULL, check_weights = TRUE)
```

## Arguments

- A:

  State matrix (`n x n`).

- B:

  Input matrix (`n x m`).

- Q:

  State weight (`n x n`, symmetric).

- R:

  Input weight (`m x m`, symmetric positive definite).

- S:

  Optional cross weight (`n x m`); defaults to zero.

- check_weights:

  If `TRUE` (default), require \\\begin{pmatrix} Q & S \\ S^\top &
  R\end{pmatrix} \succeq 0\\, the condition under which the equation
  describes a well-posed LQ problem. Set to `FALSE` to solve indefinite
  equations (for example in H-infinity synthesis).

## Value

A list of class `ctrl_riccati` with `P` (stabilizing solution), `K` (the
associated gain \\R^{-1}(B^\top P + S^\top)\\ or \\(R + B^\top P
B)^{-1}(B^\top P A + S^\top)\\), `poles` (eigenvalues of \\A - B K\\),
`residual` (relative Frobenius residual), `symmetry_error`, `rcond_u11`
and `discrete`.

## Details

Both solvers use the extended (Arnold and Laub) matrix pencil of
dimension \\2n + m\\, compress away the input block with an orthogonal
QR factorization (so `R` is never inverted explicitly), and compute an
ordered generalized real Schur decomposition with LAPACK `dgges` that
places the `n` stable eigenvalues first. The solution is \\P = U\_{21}
U\_{11}^{-1}\\, obtained with a linear solve. The discrete pencil
formulation does not require `A` to be invertible.

The stabilizing solution exists and is unique when `R` is positive
definite, \\(A, B)\\ is stabilizable and the pencil has no eigenvalues
on the stability boundary (guaranteed, for example, when \\Q - S R^{-1}
S^\top \succeq 0\\ and the corresponding pair is detectable). When the
boundary condition fails, or when \\U\_{11}\\ is singular, an error is
raised rather than an invalid solution returned. A solution whose
relative residual exceeds 1e-6, or whose closed loop is not stable,
raises an error. The returned object reports the relative residual, the
symmetry error and the closed-loop eigenvalues so that the result can be
verified.

## References

Arnold, W. F. and Laub, A. J. (1984). Generalized eigenproblem
algorithms and software for algebraic Riccati equations. *Proceedings of
the IEEE*, 72(12), 1746-1754.

## Examples

``` r
# Scalar CARE with known solution sqrt(2) - 1
ctrl_care(-1, 1, 1, 1)$P
#>           [,1]
#> [1,] 0.4142136
# Scalar DARE with known solution (1 + sqrt(5)) / 2
ctrl_dare(1, 1, 1, 1)$P
#>          [,1]
#> [1,] 1.618034
```
