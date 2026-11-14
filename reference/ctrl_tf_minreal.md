# Pole-zero cancellation for transfer functions

Remove numerator and denominator roots that coincide within `tol`
(relative to the root magnitude). Cancellations are a numerical
decision: a near cancellation of an unstable pole hides an internal
instability, so inspect the result before using it for design.

## Usage

``` r
ctrl_tf_minreal(sys, tol = sqrt(.Machine$double.eps))
```

## Arguments

- sys:

  A SISO `ctrl_tf` model.

- tol:

  Relative tolerance for matching roots.

## Value

A `ctrl_tf` model.

## Examples

``` r
G <- ctrl_tf(c(1, 1), c(1, 3, 2))  # (s + 1) / ((s + 1)(s + 2))
ctrl_tf_minreal(G)
#> Transfer function (continuous time)
#> 
#>         1
#>   ---------
#>   s +     2
```
