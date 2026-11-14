# Controllability and observability matrices

\$\$\mathcal C = \[B, AB, \ldots, A^{n-1}B\], \qquad \mathcal O = \[C;
CA; \ldots; CA^{n-1}\].\$\$ These matrices are badly scaled for all but
small systems; prefer
[`ctrl_controllability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_controllability.md)
and
[`ctrl_observability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_observability.md),
which use an orthogonal staircase reduction, for rank decisions.

## Usage

``` r
ctrl_ctrb(sys, B = NULL)

ctrl_obsv(sys, C = NULL)
```

## Arguments

- sys:

  A `ctrl_ss` model, or the matrix `A` when `B` (or `C`) is given.

- B, C:

  Optional matrices when `sys` is a matrix.

## Value

A numeric matrix.

## Examples

``` r
A <- matrix(c(0, 1, -2, -3), 2, byrow = TRUE)
ctrl_ctrb(A, c(0, 1))
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    1   -3
ctrl_obsv(A, C = c(1, 0))
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```
