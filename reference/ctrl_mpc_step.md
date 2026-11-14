# One MPC step

Solve the MPC problem at the current state and return the first input
together with the full predicted plan and solver diagnostics. Failed
solves (any OSQP status other than `"solved"`, or `"solved inaccurate"`
when `accept_inaccurate = TRUE`) are reported with `ok = FALSE` and
`u = NA`; the caller decides on a fallback (see
[`ctrl_simulate_mpc()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_simulate_mpc.md)).

## Usage

``` r
ctrl_mpc_step(
  mpc,
  x,
  reference = NULL,
  u_prev = NULL,
  disturbance = NULL,
  k = 1L,
  warm_start = NULL,
  output_offset = NULL
)
```

## Arguments

- mpc:

  A `ctrl_mpc` object.

- x:

  Current state.

- reference:

  Optional list with `x`, `u`, `y` references (vectors, matrices with
  one row per stage, or functions of the stage index).

- u_prev:

  Previously applied input (for rate terms).

- disturbance:

  Known additive disturbance \\d_k\\ in the dynamics (vector or matrix
  with `N` rows).

- k:

  Time index passed to reference functions.

- warm_start:

  Optional primal vector to warm start OSQP (for example the shifted
  previous solution).

- output_offset:

  Optional known additive output offset (length `p`), for example an
  estimated output disturbance: output bounds then apply to \\C x + D
  u + \text{offset}\\.

## Value

A list of class `ctrl_mpc_step` with `u` (input to apply), `ok`,
`status`, `x_plan` (`(N+1) x n`), `u_plan` (`N x m`), `objective`
(including constants), `slack`, `iterations`, `solve_time`,
`primal_residual`, `dual_residual`, `max_violation` (largest bound
violation of the plan, slacks excluded) and `dynamics_residual`.

## Examples

``` r
plant <- ctrl_ss(1, 1, 1, dt = 1)
mpc <- ctrl_mpc(plant, 5, 1, 1, u_min = -0.5, u_max = 0.5)
ctrl_mpc_step(mpc, x = 3)$u
#>   u1 
#> -0.5 
```
