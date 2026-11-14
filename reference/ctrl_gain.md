# Static gain model

A state-space model without states, `y = D u`.

## Usage

``` r
ctrl_gain(D, dt = NULL, input_names = NULL, output_names = NULL)
```

## Arguments

- D:

  Gain matrix (`p x m`) or scalar.

- dt:

  `NULL` or `0` for continuous time, or a positive sample interval.

- input_names:

  Optional unique input names.

- output_names:

  Optional unique output names. Defaults are `x1, x2, ...`, `u1, ...`
  and `y1, ...`.

## Value

A `ctrl_ss` object with zero states.

## Examples

``` r
ctrl_gain(diag(2))
#> State-space model (continuous time)
#> 0 states, 2 inputs, 2 outputs
#> 
#> A =
#>   (empty 0 x 0)
#> 
#> B =
#>   (empty 0 x 2)
#> 
#> C =
#>   (empty 2 x 0)
#> 
#> D =
#>    u1 u2
#> y1  1  0
#> y2  0  1
```
