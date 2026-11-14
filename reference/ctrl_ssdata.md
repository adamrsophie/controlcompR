# State-space data and dimensions

Accessors that avoid depending on the internal fields of model objects.

## Usage

``` r
ctrl_ssdata(sys)

ctrl_dims(sys)

ctrl_is_discrete(sys)

ctrl_is_continuous(sys)

ctrl_sample_time(sys)
```

## Arguments

- sys:

  A `ctrl_ss`, `ctrl_tf` or `ctrl_zpk` model (transfer functions are
  converted with
  [`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)).

## Value

`ctrl_ssdata()` returns a list with `A`, `B`, `C`, `D` and `dt`.
`ctrl_dims()` returns a named integer vector with `states`, `inputs` and
`outputs`. `ctrl_is_discrete()` and `ctrl_is_continuous()` return a
logical, and `ctrl_sample_time()` returns `0` for continuous models.

## Examples

``` r
sys <- ctrl_ss(-1, 1, 1)
ctrl_dims(sys)
#>  states  inputs outputs 
#>       1       1       1 
ctrl_ssdata(sys)$A
#>      [,1]
#> [1,]   -1
ctrl_is_discrete(ctrl_c2d(sys, 0.1))
#> [1] TRUE
```
