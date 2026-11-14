# Test for control system objects

Test for control system objects

## Usage

``` r
is_ctrl_ss(x)

is_ctrl_tf(x)

is_ctrl_system(x)
```

## Arguments

- x:

  Any object.

## Value

`TRUE` or `FALSE`.

## Examples

``` r
is_ctrl_ss(ctrl_ss(-1, 1, 1))
#> [1] TRUE
is_ctrl_tf(ctrl_tf(1, c(1, 1)))
#> [1] TRUE
```
