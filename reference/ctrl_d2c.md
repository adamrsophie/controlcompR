# Convert a discrete time model to continuous time

`"zoh"` inverts the zero-order-hold map with a matrix logarithm of
\\\begin{pmatrix} A_d & B_d \\ 0 & I\end{pmatrix}\\; it fails when
\\A_d\\ has eigenvalues on the closed negative real axis (no real
logarithm exists, the discrete model is not the sampled version of any
real continuous model). `"tustin"` inverts the bilinear map and fails
when \\A_d\\ has an eigenvalue at -1. Both are exact inverses of the
corresponding
[`ctrl_c2d()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_c2d.md)
methods when they apply.

## Usage

``` r
ctrl_d2c(sys, method = c("zoh", "tustin"))
```

## Arguments

- sys:

  A discrete time model.

- method:

  `"zoh"` or `"tustin"`.

## Value

A continuous time model.

## Examples

``` r
sysd <- ctrl_c2d(ctrl_ss(-1, 1, 1), 0.1)
ctrl_d2c(sysd)$A  # -1
#>      [,1]
#> [1,]   -1
```
