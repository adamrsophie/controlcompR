# Linearize a nonlinear model

Compute the Jacobians \\A = \partial f/\partial x\\, \\B = \partial
f/\partial u\\, \\C = \partial h/\partial x\\, \\D = \partial h/\partial
u\\ at an operating point. The result is a local model in deviation
variables \\\delta x = x - x_0\\, \\\delta u = u - u_0\\; it is exact
only to first order and only near the operating point. Jacobians are
computed by central finite differences (step \\\epsilon^{1/3}\max(1,
\|x_i\|)\\, error \\O(h^2)\\) unless analytic Jacobians are supplied.
The method used is recorded in the result.

## Usage

``` r
ctrl_linearize(
  sys,
  x0,
  u0 = numeric(sys$m),
  t = 0,
  method = c("central", "forward"),
  jacobian = NULL
)
```

## Arguments

- sys:

  A `ctrl_nlsys` model.

- x0, u0:

  Operating point.

- t:

  Time at which to linearize.

- method:

  `"central"` or `"forward"` differences.

- jacobian:

  Optional function `function(t, x, u, params)` returning
  `list(A = , B = )` (and optionally `C`, `D`) analytically.

## Value

A `ctrl_ss` model with attributes `operating_point`, `drift` and
`method`.

## Details

At an operating point that is not an equilibrium, the linearization has
an additional constant term \\f(x_0, u_0)\\ (returned as `drift`) that
the `ctrl_ss` model does not include.

## Examples

``` r
pendulum <- ctrl_nlsys(
  function(t, x, u, p) c(x[2], -9.81 * sin(x[1]) - 0.5 * x[2] + u),
  output = function(t, x, u, p) x[1], nstates = 2, ninputs = 1
)
lin <- ctrl_linearize(pendulum, x0 = c(pi, 0), u0 = 0)
ctrl_poles(lin)  # upright position: unstable
#> [1]  2.892053+0i -3.392053+0i
```
