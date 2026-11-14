# Nonlinear state-space models

Define a nonlinear model \\\dot x = f(t, x, u, \theta)\\ (continuous
time, `dt = NULL`) or \\x\_{k+1} = f(k, x_k, u_k, \theta)\\ (discrete
time, `dt > 0`) with output \\y = h(t, x, u, \theta)\\.

## Usage

``` r
ctrl_nlsys(
  update,
  output = NULL,
  nstates,
  ninputs,
  noutputs = NULL,
  dt = NULL,
  params = list(),
  state_names = NULL,
  input_names = NULL,
  output_names = NULL
)
```

## Arguments

- update:

  Function `function(t, x, u, params)` returning the state derivative
  (continuous) or the next state (discrete).

- output:

  Function `function(t, x, u, params)` returning the output; `NULL`
  outputs the full state.

- nstates, ninputs:

  Dimensions.

- noutputs:

  Number of outputs (defaults to `nstates` without `output`; otherwise
  evaluated once).

- dt:

  `NULL` for continuous time or a positive sample interval.

- params:

  List of parameters passed to `update` and `output`.

- state_names, input_names, output_names:

  Optional names.

## Value

An object of class `ctrl_nlsys`.

## See also

[`ctrl_nl_simulate()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nl_simulate.md),
[`ctrl_linearize()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_linearize.md),
[`ctrl_equilibrium()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_equilibrium.md)

## Examples

``` r
pendulum <- ctrl_nlsys(
  function(t, x, u, p) c(x[2], -p$g / p$l * sin(x[1]) - p$c * x[2] + u),
  output = function(t, x, u, p) x[1],
  nstates = 2, ninputs = 1, params = list(g = 9.81, l = 1, c = 0.5)
)
pendulum
#> Nonlinear model (continuous time): 2 states, 1 inputs, 1 outputs
```
