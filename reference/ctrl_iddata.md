# Identification data

Store sampled input-output data for system identification.

## Usage

``` r
ctrl_iddata(y, u, dt = 1, detrend = FALSE)
```

## Arguments

- y:

  Output data (vector or `T x p` matrix).

- u:

  Input data (vector or `T x m` matrix).

- dt:

  Sample interval.

- detrend:

  Remove the sample means of `y` and `u` (stored in `offsets` so
  predictions can be mapped back), or a list with `y` and `u` offsets to
  subtract, typically `train$offsets`. Models estimated on detrended
  data describe deviations from the offsets, so validation data must be
  shifted by the training offsets, not by their own means.

## Value

An object of class `ctrl_iddata`.

## Examples

``` r
d <- ctrl_iddata(y = sin(1:50 / 5), u = cos(1:50 / 5), dt = 0.1)
d
#> Identification data: 50 samples, 1 input(s), 1 output(s), dt = 0.1
```
