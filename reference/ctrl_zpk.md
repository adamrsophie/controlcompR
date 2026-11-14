# Zero-pole-gain models

A SISO model \\G(s) = k \prod (s - z_i) / \prod (s - p_i)\\. Complex
zeros and poles must appear in conjugate pairs so that the model has
real coefficients.

## Usage

``` r
ctrl_zpk(zeros = numeric(), poles = numeric(), gain = 1, dt = NULL)
```

## Arguments

- zeros:

  Vector of zeros (possibly complex).

- poles:

  Vector of poles (possibly complex).

- gain:

  Real gain `k`.

- dt:

  `NULL` or `0` for continuous time, or a positive sample interval.

## Value

An object of class `c("ctrl_zpk", "ctrl_system")`.

## Examples

``` r
H <- ctrl_zpk(-1, c(-2, -3 + 1i, -3 - 1i), 4)
H
#> Zero-pole-gain model (continuous time)
#>   gain: 4
#>   zeros: -1+0i
#>   poles: -2+0i, -3+1i, -3-1i
ctrl_zpk2tf(H)
#> Transfer function (continuous time)
#> 
#>                4 s +     4
#>   ---------------------------------
#>   s^3 +     8 s^2 +    22 s +    20
```
