# Simulate a nonlinear model

Continuous time models are integrated with `deSolve` (default method
`"lsoda"`, an adaptive stiff/nonstiff solver); discrete models are
iterated. The input can be

## Usage

``` r
ctrl_nl_simulate(
  sys,
  x0,
  times,
  input = NULL,
  controller = NULL,
  method = "lsoda",
  rtol = 1e-08,
  atol = 1e-10
)
```

## Arguments

- sys:

  A `ctrl_nlsys` model.

- x0:

  Initial state.

- times:

  Time points (for discrete models, sample indices are `times / dt`).

- input:

  Input (see above) or `NULL` for zero input.

- controller:

  Optional feedback law `function(t, x)`.

- method:

  `deSolve` integration method.

- rtol, atol:

  Integration tolerances.

## Value

An object of class `ctrl_nl_response` with `time`, `x`, `u` and `y`
matrices.

## Details

- a vector or matrix with one row per time point, held constant between
  points (zero-order hold, integrated interval by interval);

- a function `input(t)` evaluated continuously by the integrator;

- a feedback law `controller(t, x)` evaluated at the time points and
  held in between, which simulates a sampled-data controller.

## Examples

``` r
pendulum <- ctrl_nlsys(
  function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
  nstates = 2, ninputs = 1
)
sim <- ctrl_nl_simulate(pendulum, x0 = c(2, 0), times = seq(0, 10, 0.05))
plot(sim)
```
