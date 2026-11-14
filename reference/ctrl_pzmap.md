# Pole-zero map data

Pole-zero map data

## Usage

``` r
ctrl_pzmap(sys)
```

## Arguments

- sys:

  A model.

## Value

An object of class `ctrl_pzmap` with `poles`, `zeros` and `dt`.

## Examples

``` r
ctrl_pzmap(ctrl_tf(c(1, 2), c(1, 2, 5)))
#> Poles: -1+2i, -1-2i
#> Zeros: -2+0i
```
