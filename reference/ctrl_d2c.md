# Convert a discrete time model to continuous time

`"zoh"` inverts the zero-order-hold map with a matrix logarithm of
\\\begin{pmatrix} A_d & B_d \\ 0 & I\end{pmatrix}\\; it fails when
\\A_d\\ has eigenvalues on the closed negative real axis. Such matrices
may still have real logarithms (for example from pairs of Jordan
blocks), but the principal logarithm is not real and other branches are
not supported. `"tustin"` inverts the bilinear map and fails when
\\A_d\\ has an eigenvalue at -1. Both are exact inverses of the
corresponding
[`ctrl_c2d()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_c2d.md)
methods when they apply.

## Usage

``` r
ctrl_d2c(sys, method = c("zoh", "tustin"), prewarp = NULL)
```

## Arguments

- sys:

  A discrete time model.

- method:

  `"zoh"` or `"tustin"`.

- prewarp:

  Prewarping frequency used in
  [`ctrl_c2d()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_c2d.md)
  (`"tustin"` only).

## Value

A continuous time model.

## Examples

``` r
sysd <- ctrl_c2d(ctrl_ss(-1, 1, 1), 0.1)
ctrl_d2c(sysd)$A  # -1
#>      [,1]
#> [1,]   -1
```
