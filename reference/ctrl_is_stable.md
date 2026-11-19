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
  by the norm of the balanced `A`. For a model whose time scales differ
  by about 1e13 or more, the slowest stable poles can fall inside this
  band and are classified as marginal; pass a smaller `tol` when such
  poles are known to be accurate.

## Value

`TRUE` when the model is asymptotically stable.

## Examples

``` r
ctrl_is_stable(ctrl_ss(-1, 1, 1))
#> [1] TRUE
```
