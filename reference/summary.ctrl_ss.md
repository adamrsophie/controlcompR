# Summaries of control models

Summaries of control models

## Usage

``` r
# S3 method for class 'ctrl_ss'
summary(object, ...)
```

## Arguments

- object:

  A model.

- ...:

  Unused.

## Value

A `ctrl_summary` list with dimensions, time domain, poles, stability,
controllability and observability ranks.

## Examples

``` r
summary(ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#> Model summary (continuous time)
#>   states: 2, inputs: 1, outputs: 1
#>   stability: asymptotically stable
#>   poles: -1+0i, -2+0i
#>   controllability rank: 2 of 2
#>   observability rank: 2 of 2
#>   DC gain:
#>     u1
#> y1 0.5
```
