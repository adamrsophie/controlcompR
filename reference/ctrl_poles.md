# Poles of a model

For a state-space model the poles are the eigenvalues of `A`. For a
transfer function they are the roots of the denominator polynomials
(after the coefficient normalization done by
[`ctrl_tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf.md));
pole-zero cancellations are not removed.

## Usage

``` r
ctrl_poles(sys)
```

## Arguments

- sys:

  A `ctrl_ss`, `ctrl_tf` or `ctrl_zpk` model.

## Value

A complex vector of poles sorted by decreasing real part (continuous
time) or decreasing modulus (discrete time).

## Examples

``` r
ctrl_poles(ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#> [1] -1+0i -2+0i
ctrl_poles(ctrl_tf(1, c(1, 3, 2)))
#> [1] -1+0i -2+0i
```
