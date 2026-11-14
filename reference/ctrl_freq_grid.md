# Frequency grid for a model

A logarithmic grid spanning two decades beyond the smallest and largest
nonzero pole and zero magnitudes, refined near lightly damped poles. For
discrete models the grid ends at the Nyquist frequency.

## Usage

``` r
ctrl_freq_grid(sys, n = 400L)
```

## Arguments

- sys:

  A model.

- n:

  Number of points in the base grid.

## Value

A sorted vector of positive frequencies.

## Examples

``` r
range(ctrl_freq_grid(ctrl_tf(1, c(1, 0.2, 100))))
#> [1]   0.1 100.0
```
