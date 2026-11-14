# Sylvester equation

Solve \\A X + X B = C\\ with the Bartels-Stewart algorithm. A unique
solution exists when `A` and `-B` have no common eigenvalues.

## Usage

``` r
ctrl_sylvester(A, B, C, tol = 1e-12)
```

## Arguments

- A:

  Square `n x n` matrix.

- B:

  Square `m x m` matrix.

- C:

  `n x m` matrix.

- tol:

  Relative uniqueness tolerance.

## Value

The `n x m` solution with a `residual` attribute.

## Examples

``` r
X <- ctrl_sylvester(diag(c(1, 2)), matrix(3), matrix(1, 2, 1))
X
#>      [,1]
#> [1,] 0.25
#> [2,] 0.20
#> attr(,"residual")
#> [1] 0
```
