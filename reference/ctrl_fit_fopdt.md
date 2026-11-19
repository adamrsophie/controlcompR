# Fit a first-order plus dead-time model to a step response

Estimate \\K e^{-\theta s}/(\tau s + 1)\\ from a step response (the gain
is divided by the step amplitude) with the two-point method (28.3
percent and 63.2 percent of the final change). The response must be
monotone enough for the crossing times to be meaningful.

## Usage

``` r
ctrl_fit_fopdt(response)
```

## Arguments

- response:

  A `ctrl_response` from
  [`ctrl_step()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md),
  or a model (its step response is computed).

## Value

A named vector `c(K, tau, theta)`.

## Examples

``` r
ctrl_fit_fopdt(ctrl_tf(1, c(1, 3, 3, 1)))
#>        K      tau    theta 
#>  1.00000 86.53948 17.93607 
```
