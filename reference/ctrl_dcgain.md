# Steady-state (DC) gain

The gain at zero frequency: \\D - C A^{-1} B\\ in continuous time and
\\D + C (I - A)^{-1} B\\ in discrete time. An error is raised when the
input-output map has a pole at `s = 0` (or `z = 1`), because the gain is
then infinite or undefined; poles there that are uncontrollable or
unobservable are removed first.

## Usage

``` r
ctrl_dcgain(sys)
```

## Arguments

- sys:

  A model.

## Value

A `p x m` matrix.

## Examples

``` r
ctrl_dcgain(ctrl_tf(2, c(1, 3, 2)))
#>    u1
#> y1  1
```
