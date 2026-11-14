# Compare controllers on a common scenario

Run several closed-loop simulation functions on the same plant and
scenario and collect their performance metrics in one table. Each
element of `controllers` is a function of no arguments that returns a
list (or data frame) with `time`, `y`, and optionally `u` and `r`; this
keeps the comparison independent of the controller type (PID, LQR, MPC
or user code).

## Usage

``` r
ctrl_compare_controllers(controllers, reference = 0, u_min = -Inf, u_max = Inf)
```

## Arguments

- controllers:

  Named list of simulation functions.

- reference:

  Reference used when a result has no `r`.

- u_min, u_max:

  Input limits for violation reporting.

## Value

A data frame with one row per controller.

## Examples

``` r
plant <- ctrl_tf(1, c(1, 1, 0))
sims <- list(
  pid_slow = function() ctrl_simulate_pid(plant, ctrl_pid(1, 0.1, 1,
    derivative_filter = 10, dt = 0.05), steps = 300),
  pid_fast = function() ctrl_simulate_pid(plant, ctrl_pid(4, 0.5, 2,
    derivative_filter = 10, dt = 0.05), steps = 300)
)
ctrl_compare_controllers(sims, reference = 1)
#>   controller      IAE      ISE      ITAE      RMSE control_effort max_abs_u
#> 1   pid_slow 2.813899 1.249996 12.465280 0.2915594      0.5989531         1
#> 2   pid_fast 1.299635 0.525195  4.139549 0.1915204      3.8672182         4
#>   max_violation
#> 1             0
#> 2             0
```
