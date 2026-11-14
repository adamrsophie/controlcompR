# Test whether a transfer function is proper

Test whether a transfer function is proper

## Usage

``` r
ctrl_is_proper(sys, strict = FALSE)
```

## Arguments

- sys:

  A `ctrl_tf` model.

- strict:

  If `TRUE`, test strict properness.

## Value

`TRUE` or `FALSE`.

## Examples

``` r
ctrl_is_proper(ctrl_tf(c(1, 0), c(1, 1)))
#> [1] TRUE
ctrl_is_proper(ctrl_tf(c(1, 0), c(1, 1)), strict = TRUE)
#> [1] FALSE
```
