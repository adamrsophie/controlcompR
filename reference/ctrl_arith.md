# Arithmetic on models

`+` is a parallel connection, `-` subtracts (or negates), `*` is a
series connection with MATLAB semantics (`G2 * G1` applies `G1` first),
`/` divides by a scalar or (SISO transfer functions) multiplies by the
inverse, and `^` raises a model to a nonnegative integer power (negative
powers for SISO transfer functions).

## Usage

``` r
# S3 method for class 'ctrl_system'
Ops(e1, e2)
```

## Arguments

- e1, e2:

  Models or numbers.

## Value

A model.

## Examples

``` r
s <- ctrl_tf("s")
G <- 1 / (s + 1)
2 * G + G^2
#> Transfer function (continuous time)
#> 
#>          2 s^2 +     5 s +     3
#>   ---------------------------------
#>   s^3 +     3 s^2 +     3 s +     1
```
