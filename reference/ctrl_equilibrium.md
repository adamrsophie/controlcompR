# Equilibrium points

Solve \\f(x, u) = 0\\ (continuous) or \\f(x, u) = x\\ (discrete) with a
damped Newton method using finite-difference Jacobians. By default the
inputs are fixed at `u0` and all states are unknown; with `free_inputs`,
the listed input components are also unknowns, and with `fixed_states`,
the listed states are held at their values in `x0` (for example to find
the input that holds a given operating point). When there are more
unknowns than equations, minimum-norm Newton steps are used and the
solution depends on the initial guess. Local stability is assessed from
the linearization: an equilibrium with a linearization on the stability
boundary is reported as `"inconclusive"`, since nonlinear terms decide
stability there.

## Usage

``` r
ctrl_equilibrium(
  sys,
  x0,
  u0 = numeric(sys$m),
  free_inputs = integer(),
  fixed_states = integer(),
  tol = 1e-10,
  maxit = 100L
)
```

## Arguments

- sys:

  A `ctrl_nlsys` model.

- x0:

  Initial guess for the state.

- u0:

  Input (or initial guess for the free inputs).

- free_inputs:

  Indices of inputs to solve for.

- fixed_states:

  Indices of states held at `x0`.

- tol:

  Residual tolerance.

- maxit:

  Maximum Newton iterations.

## Value

A list with `x`, `u`, `residual`, `converged`, `iterations`,
`eigenvalues` (of the linearization) and `stability`.

## Examples

``` r
pendulum <- ctrl_nlsys(
  function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
  nstates = 2, ninputs = 1
)
ctrl_equilibrium(pendulum, x0 = c(3, 0))$x
#> [1]  3.141593e+00 -2.587262e-25
```
