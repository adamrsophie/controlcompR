# Step-response characteristics

Compute the usual transient characteristics for each output and input
case of a step response. Definitions (with \\y_0\\ the value before the
step and \\y\_\infty\\ the final value; when the model is known, a
direct feedthrough jump at \\t = 0\\ is not counted in \\y_0\\):

## Usage

``` r
ctrl_step_info(response, sys = NULL, settle = 0.02, rise = c(0.1, 0.9))
```

## Arguments

- response:

  A `ctrl_response` (typically from
  [`ctrl_step()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md))
  or a model, whose step response is then computed.

- sys:

  Optional model for the exact final value.

- settle:

  Relative settling band.

- rise:

  Rise-time fractions.

## Value

A data frame with one row per output and case.

## Details

- rise time: time to go from `rise[1]` to `rise[2]` (default 10 to 90
  percent) of \\y\_\infty - y_0\\, with linear interpolation;

- settling time: the last time the response is outside the band \\\|y -
  y\_\infty\| \le\\ `settle` \\\cdot \|y\_\infty - y_0\|\\ (default 2
  percent); `NA` if it never enters the band for good on the simulated
  horizon;

- overshoot and undershoot: percent of \\\|y\_\infty - y_0\|\\;

- peak and peak time: largest \\\|y - y_0\|\\.

When `sys` (the model) is supplied the final value is its DC gain times
the step amplitude (exact, for stable models). Otherwise the last sample
is used and `settled` reports whether the response varies by less than
the band over the last 5 percent of the horizon. Characteristics depend
on these conventions, so state them when reporting numbers.

## Examples

``` r
G <- ctrl_tf(4, c(1, 1.2, 4))
ctrl_step_info(G)
#>   output    case rise_time settling_time overshoot undershoot     peak
#> 1     y1 step u1 0.6607187      5.623333  37.23241          0 1.372324
#>   peak_time final_value settled
#> 1     1.645           1    TRUE
```
