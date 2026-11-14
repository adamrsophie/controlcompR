# Bandwidth and resonant peak

`ctrl_bandwidth()` returns the first frequency at which the gain of a
SISO model drops `drop_db` (default 3 dB) below its DC gain; it is `NA`
when the gain never drops that far on the analyzed range and an error
when the DC gain is infinite. `ctrl_peak()` returns the peak gain over
frequency, which equals the H-infinity norm for stable models.

## Usage

``` r
ctrl_bandwidth(sys, drop_db = 3)

ctrl_peak(sys)
```

## Arguments

- sys:

  A SISO model.

- drop_db:

  Gain drop in dB.

## Value

`ctrl_bandwidth()`: a number. `ctrl_peak()`: a list with `gain` and
`frequency`.

## Examples

``` r
G <- ctrl_tf(1, c(1, 0.5, 1))
ctrl_bandwidth(G)
#> [1] 1.483908
ctrl_peak(G)
#> $gain
#> [1] 2.065591
#> 
#> $frequency
#> [1] 0.9354142
#> 
```
