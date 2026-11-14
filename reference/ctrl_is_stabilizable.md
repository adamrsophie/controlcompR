# Stabilizability and detectability tests

Stabilizability and detectability tests

## Usage

``` r
ctrl_is_stabilizable(sys, tol = NULL)

ctrl_is_detectable(sys, tol = NULL)

ctrl_is_controllable(sys, tol = NULL)

ctrl_is_observable(sys, tol = NULL)
```

## Arguments

- sys:

  A `ctrl_ss` model (or the matrix `A`).

- tol:

  Rank tolerance.

## Value

`TRUE` or `FALSE`.

## Examples

``` r
ctrl_is_stabilizable(ctrl_ss(diag(c(-1, 2)), c(0, 1), diag(2)))
#> [1] TRUE
ctrl_is_detectable(ctrl_ss(diag(c(-1, 2)), c(0, 1), c(1, 0)))
#> [1] FALSE
```
