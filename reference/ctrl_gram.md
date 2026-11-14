# Controllability and observability Gramians

For an asymptotically stable model, the controllability Gramian solves
\\A W_c + W_c A^\top + B B^\top = 0\\ (continuous) or \\A W_c A^\top -
W_c + B B^\top = 0\\ (discrete), and the observability Gramian solves
the dual equations with \\A^\top\\ and \\C^\top C\\.

## Usage

``` r
ctrl_gram(sys, type = c("c", "o"))
```

## Arguments

- sys:

  A `ctrl_ss` model (asymptotically stable).

- type:

  `"c"` (controllability) or `"o"` (observability).

## Value

A symmetric positive semidefinite matrix.

## Examples

``` r
sys <- ctrl_ss(matrix(c(-1, 0, 1, -2), 2), c(0, 1), c(1, 0))
ctrl_gram(sys, "c")
#>            [,1]       [,2]
#> [1,] 0.08333333 0.08333333
#> [2,] 0.08333333 0.25000000
#> attr(,"residual")
#> [1] 3.925231e-17
#> attr(,"separation")
#> [1] 1
```
