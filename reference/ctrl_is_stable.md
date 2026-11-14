# Test asymptotic stability

Test asymptotic stability

## Usage

``` r
ctrl_is_stable(sys, tol = NULL)
```

## Arguments

- sys:

  A model.

- tol:

  Boundary tolerance. Defaults to a multiple of machine precision scaled
  by the norm of `A`.

## Value

`TRUE` when the model is asymptotically stable.

## Examples

``` r
ctrl_is_stable(ctrl_ss(-1, 1, 1))
#> [1] TRUE
```
