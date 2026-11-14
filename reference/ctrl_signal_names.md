# Signal names of a model

Signal names of a model

## Usage

``` r
ctrl_signal_names(sys)
```

## Arguments

- sys:

  A model.

## Value

A list with `states`, `inputs` and `outputs` character vectors.

## Examples

``` r
ctrl_signal_names(ctrl_ss(-1, 1, 1, input_names = "force"))
#> $states
#> [1] "x1"
#> 
#> $inputs
#> [1] "force"
#> 
#> $outputs
#> [1] "y1"
#> 
```
