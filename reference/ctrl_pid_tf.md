# Linear transfer functions of a PID controller

Return the linear part of the controller (ignoring saturation) as the
transfer function from the error \\e = r - y\\ to `u` (this one-input
form corresponds to setpoint weights `b = c = 1`), \\C(s) = k_p +
k_i/s + k_d s/(T_f s + 1)\\, or, with `two_degree = TRUE`, as a
two-input model from `(r, y)` to `u` that includes the setpoint weights.
When the controller has a sample interval, the discrete equivalent of
the implemented difference equations is returned.

## Usage

``` r
ctrl_pid_tf(pid, two_degree = FALSE)
```

## Arguments

- pid:

  A `ctrl_pid` object.

- two_degree:

  Return the two-input form.

## Value

A `ctrl_tf` model.

## Examples

``` r
ctrl_pid_tf(ctrl_pid(kp = 1, ki = 2))
#> Transfer function (continuous time)
#> 
#>   s +     2
#>   ---------
#>       s
```
