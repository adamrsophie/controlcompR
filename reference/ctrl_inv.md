# Inverse of a model

SISO transfer functions are inverted by swapping numerator and
denominator (the result may be improper). State-space models need a
square, invertible feedthrough `D`; the inverse is \\(A - B D^{-1} C, B
D^{-1}, -D^{-1} C, D^{-1})\\.

## Usage

``` r
ctrl_inv(sys)
```

## Arguments

- sys:

  A model.

## Value

The inverse model.

## Examples

``` r
ctrl_inv(ctrl_tf(c(1, 2), c(1, 3)))
#> Transfer function (continuous time)
#> 
#>   s +     3
#>   ---------
#>   s +     2
```
