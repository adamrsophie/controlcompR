# Pole placement

Compute a state-feedback gain `K` such that the eigenvalues of \\A - B
K\\ are the requested `poles`.

## Usage

``` r
ctrl_place(sys, B = NULL, poles, maxiter = 50L)
```

## Arguments

- sys:

  A `ctrl_ss` model or the matrix `A`.

- B:

  Input matrix when `sys` is a matrix.

- poles:

  Desired closed-loop poles (length `n`).

- maxiter:

  Maximum number of KNV sweeps.

## Value

The `m x n` gain matrix `K` with attributes `poles` (achieved), `cond_X`
and `max_error` (largest distance between requested and achieved poles).

## Details

The method is the robust eigenstructure assignment of Kautsky, Nichols
and Van Dooren (method 0). With \\B = U_0 Z\\ (QR), each closed-loop
eigenvector must lie in \\S_i = \mathrm{null}(U_1^\top (A - \lambda_i
I))\\; the eigenvectors are chosen iteratively to make the eigenvector
matrix `X` as well conditioned as possible, and \\K = Z^{-1} U_0^\top
(A - X \Lambda X^{-1})\\. For single-input systems the eigenvectors are
unique, so the result equals Ackermann's formula, but it is computed
without the controllability matrix.

Requirements: `(A, B)` controllable, `B` of full column rank, complex
poles in conjugate pairs, and no pole repeated more than `rank(B)`
times. The achieved closed-loop poles are verified; `cond_X` (condition
number of the eigenvector matrix) measures how sensitive the placement
is to perturbations; large values mean the placed poles are fragile.

## References

Kautsky, J., Nichols, N. K. and Van Dooren, P. (1985). Robust pole
assignment in linear state feedback. *International Journal of Control*,
41(5), 1129-1155.

## Examples

``` r
plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
K <- ctrl_place(plant, poles = c(-4, -5))
eigen(plant$A - plant$B %*% K)$values
#> [1] -5 -4
```
