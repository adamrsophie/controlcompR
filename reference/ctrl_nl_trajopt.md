# Nonlinear trajectory optimization

Compute an open-loop input sequence \\u_0, \ldots, u\_{N-1}\\ for a
discrete time `ctrl_nlsys` that minimizes \$\$\sum\_{k=0}^{N-1} \\x_k -
x^r_k\\\_Q^2 + \\u_k - u^r_k\\\_R^2 + \\x_N - x^r_N\\\_P^2 + \rho \sum_k
\\\max(0, x_k - x\_{\max}, x\_{\min} - x_k)\\^2\$\$ from a given initial
state, by single shooting with L-BFGS-B and the adjoint gradient of
[`ctrl_nmpc()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nmpc.md).
Input bounds are hard; state bounds are penalized.

## Usage

``` r
ctrl_nl_trajopt(
  sys,
  x0,
  horizon,
  Q,
  R,
  P = NULL,
  x_ref = NULL,
  u_ref = NULL,
  x_final = NULL,
  terminal_weight = 100,
  terminal_tol = 1e-06,
  max_outer = 30L,
  u_min = -Inf,
  u_max = Inf,
  x_min = -Inf,
  x_max = Inf,
  penalty = 10000,
  maxit = 1000L,
  u_init = NULL
)
```

## Arguments

- sys:

  Discrete time `ctrl_nlsys`.

- x0:

  Initial state.

- horizon:

  Prediction horizon.

- Q, R, P:

  Weights (`P` defaults to `Q`).

- x_ref, u_ref:

  References (see
  [`ctrl_nmpc_step()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nmpc.md)).

- x_final:

  Optional target state at the final sample.

- terminal_weight:

  Penalty parameter \\\rho\\ for `x_final`.

- terminal_tol:

  Required accuracy of the terminal state.

- max_outer:

  Maximum number of multiplier updates.

- u_min, u_max:

  Input bounds.

- x_min, x_max:

  State bounds (soft, penalized).

- penalty:

  Penalty weight for state bound violations.

- maxit:

  Maximum optimizer iterations per step.

- u_init:

  Initial guess (`N x m` matrix).

## Value

An object of class `ctrl_trajopt` with `time`, `x` (`N + 1` rows), `u`
(`N` rows), `objective` (the cost without terminal constraint terms),
`converged`, `convergence` (optimizer code of the last subproblem),
`message`, `outer_iterations`, `multiplier` and `terminal_error`
(distance to `x_final`, `NA` without one).

## Details

**Terminal constraint.** A target `x_final` is imposed as the equality
\\x_N = x_f\\ by an augmented Lagrangian method: each outer iteration
minimizes the cost plus \\\nu^\top (x_N - x_f) + \frac{\rho}{2}\\x_N -
x_f\\^2\\ and then updates \\\nu \leftarrow \nu + \rho (x_N - x_f)\\
(multiplying \\\rho\\ by 10 when the violation shrinks by less than a
factor 4), until \\\\x_N - x_f\\ \le\\ `terminal_tol`. A moderate
\\\rho\\ (`terminal_weight`) keeps each subproblem well conditioned,
unlike a single large penalty. `P` is ignored when `x_final` is given.

The problem is nonconvex and the result is a local minimum that depends
on `u_init`. Failures are reported in `converged` and `message`, never
hidden.

## Examples

``` r
f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] + 0.1 * (-sin(x[1]) + u))
sys <- ctrl_nlsys(f, nstates = 2, ninputs = 1, dt = 0.1)
tr <- ctrl_nl_trajopt(sys, x0 = c(0, 0), horizon = 40, Q = 0, R = 0.01,
                      x_final = c(pi / 2, 0), u_min = -2, u_max = 2)
tr$terminal_error
#> [1] 9.058934e-07
```
