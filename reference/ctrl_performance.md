# Tracking and control performance metrics

Integral criteria for a tracking or regulation experiment on a time grid
(trapezoidal rule in continuous time, sums times `dt` for sampled data):
IAE \\\int \|e\|\\, ISE \\\int e^2\\, ITAE \\\int t\|e\|\\, RMSE
\\\sqrt{\mathrm{mean}(e^2)}\\, control effort \\\int u^2\\, maximum
\\\|u\|\\ and, when limits are given, the largest constraint violation.

## Usage

``` r
ctrl_performance(
  time,
  y,
  r = 0,
  u = NULL,
  u_min = -Inf,
  u_max = Inf,
  discrete = FALSE
)
```

## Arguments

- time:

  Time vector.

- y:

  Output matrix or vector.

- r:

  Reference (same shape as `y`, or a scalar).

- u:

  Optional control signal.

- u_min, u_max:

  Optional input limits for violation reporting.

- discrete:

  Use rectangle sums instead of the trapezoidal rule.

## Value

A one-row data frame.

## Examples

``` r
t <- seq(0, 10, by = 0.01)
ctrl_performance(t, 1 - exp(-t), 1)
#>         IAE       ISE      ITAE      RMSE
#> 1 0.9999629 0.5000167 0.9994923 0.2246135
```
