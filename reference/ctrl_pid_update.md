# One step of a discrete PID controller

Advance the discrete PID implementation by one sample. The state holds
the integrator, the filtered derivative, the previous derivative input
and the mode.

## Usage

``` r
ctrl_pid_update(pid, state = NULL, r, y, manual = NULL)
```

## Arguments

- pid:

  A `ctrl_pid` object with a positive `dt`.

- state:

  Controller state from the previous call (`NULL` to initialize).

- r, y:

  Reference and measurement at this sample.

- manual:

  Optional manual output; when given, the controller is in manual mode
  and tracks this value for bumpless transfer.

## Value

A list with `u` (saturated output), `u_unsat`, the parts `P`, `I`, `D`,
`saturated` and the updated `state`.

## Examples

``` r
pid <- ctrl_pid(kp = 1, ki = 1, dt = 0.1)
st <- NULL
for (k in 1:3) {
  out <- ctrl_pid_update(pid, st, r = 1, y = 0)
  st <- out$state
}
out$u
#> [1] 1.2
```
