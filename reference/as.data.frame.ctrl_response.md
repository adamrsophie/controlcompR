# Long data frame of a time response

Long data frame of a time response

## Usage

``` r
# S3 method for class 'ctrl_response'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `ctrl_response`.

- row.names, optional:

  Unused.

- ...:

  Unused.

## Value

A data frame with columns `time`, `case`, `kind` (`"output"`, `"state"`
or `"input"`), `signal` and `value`.

## Examples

``` r
as.data.frame(ctrl_step(ctrl_tf(1, c(1, 1)), times = 0:3))
#>    time    case   kind signal     value
#> 1     0 step u1 output     y1 0.0000000
#> 2     1 step u1 output     y1 0.6321206
#> 3     2 step u1 output     y1 0.8646647
#> 4     3 step u1 output     y1 0.9502129
#> 5     0 step u1  state     x1 0.0000000
#> 6     1 step u1  state     x1 0.6321206
#> 7     2 step u1  state     x1 0.8646647
#> 8     3 step u1  state     x1 0.9502129
#> 9     0 step u1  input     u1 1.0000000
#> 10    1 step u1  input     u1 1.0000000
#> 11    2 step u1  input     u1 1.0000000
#> 12    3 step u1  input     u1 1.0000000
```
