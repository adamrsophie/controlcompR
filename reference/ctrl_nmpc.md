# Nonlinear model predictive control (experimental)

Single-shooting nonlinear MPC for a discrete time `ctrl_nlsys`: at each
sample the input sequence minimizing \$\$\sum\_{k=0}^{N-1} \\x_k -
x^r\\\_Q^2 + \\u_k - u^r\\\_R^2 + \\x_N - x^r\\\_P^2 + \rho \sum_k
\\\max(0, x_k - x\_{\max}, x\_{\min} - x_k)\\^2\$\$ over the box
\\u\_{\min} \le u_k \le u\_{\max}\\ is computed with the bounded
quasi-Newton method L-BFGS-B
([`stats::optim`](https://rdrr.io/r/stats/optim.html)). The gradient is
computed by a backward (adjoint) pass through the prediction, with the
stage Jacobians of the model approximated by central differences, so it
costs \\O(N (n + m))\\ model evaluations. Input bounds are hard; state
bounds are enforced only through the quadratic penalty \\\rho\\ and can
be violated.

## Usage

``` r
ctrl_nmpc(
  sys,
  horizon,
  Q,
  R,
  P = NULL,
  u_min = -Inf,
  u_max = Inf,
  x_min = -Inf,
  x_max = Inf,
  penalty = 10000,
  maxit = 200L
)

ctrl_nmpc_step(nmpc, x, x_ref = NULL, u_ref = NULL, u_init = NULL, k = 0L)

ctrl_simulate_nmpc(
  nmpc,
  x0,
  steps,
  x_ref = NULL,
  u_ref = NULL,
  plant = NULL,
  on_failure = c("stop", "hold")
)
```

## Arguments

- sys:

  Discrete time `ctrl_nlsys`.

- horizon:

  Prediction horizon.

- Q, R, P:

  Weights (`P` defaults to `Q`).

- u_min, u_max:

  Input bounds.

- x_min, x_max:

  State bounds (soft, penalized).

- penalty:

  Penalty weight for state bound violations.

- maxit:

  Maximum optimizer iterations per step.

- nmpc:

  A `ctrl_nmpc` object.

- x:

  Current state.

- x_ref, u_ref:

  State and input references: vectors (constant) or matrices whose row
  `i` is the reference at sample `i - 1` (the last row is repeated
  beyond the end), or functions of the 1-based sample index.

- u_init:

  Initial guess (`N x m` matrix).

- k:

  Time index.

- x0:

  Initial state.

- steps:

  Number of closed-loop samples.

- plant:

  Simulation model (default: the prediction model).

- on_failure:

  `"stop"` or `"hold"` (apply the previous input).

## Value

An object of class `ctrl_nmpc`.

## Details

**Limitations.** The problem is nonconvex: the solver returns a local
minimum that depends on the initial guess (the shifted previous solution
is used), and no stability or recursive feasibility guarantee is
provided. Each step reports the optimizer convergence code; nonzero
codes are treated as failures.

## Examples

``` r
f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] + 0.1 * (-sin(x[1]) + u))
sys <- ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1)
nmpc <- ctrl_nmpc(sys, horizon = 15, Q = diag(c(10, 1)), R = 0.1,
                  u_min = -2, u_max = 2)
sim <- ctrl_simulate_nmpc(nmpc, x0 = c(1, 0), steps = 30)
tail(sim$x, 1)
#>                  x1          x2
#> [31,] -0.0008910959 0.003059755
```
