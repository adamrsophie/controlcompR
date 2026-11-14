# Hankel singular values

The square roots of the eigenvalues of \\W_c W_o\\ for an asymptotically
stable model, computed from Gramian factors as the singular values of
\\L_o^\top L_c\\.

## Usage

``` r
ctrl_hsvd(sys)
```

## Arguments

- sys:

  An asymptotically stable model.

## Value

A decreasing numeric vector.

## Examples

``` r
ctrl_hsvd(ctrl_tf(1, c(1, 3, 2)))
#> [1] 0.29679607 0.04679607
```
