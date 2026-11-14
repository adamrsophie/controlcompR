# Discretize a continuous nonlinear model

Build a discrete time `ctrl_nlsys` whose update integrates the
continuous model over one sample interval with the input held constant
(`deSolve`, tight tolerances, or a fixed number of classical Runge-Kutta
steps with `method = "rk4"`).

## Usage

``` r
ctrl_nl_discretize(sys, dt, method = c("ode", "rk4"), substeps = 4L)
```

## Arguments

- sys:

  Continuous `ctrl_nlsys`.

- dt:

  Sample interval.

- method:

  `"ode"` (deSolve lsoda) or `"rk4"`.

- substeps:

  RK4 steps per sample.

## Value

A discrete `ctrl_nlsys`.

## Examples

``` r
pendulum <- ctrl_nlsys(
  function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
  nstates = 2, ninputs = 1
)
pd <- ctrl_nl_discretize(pendulum, 0.05, method = "rk4")
```
