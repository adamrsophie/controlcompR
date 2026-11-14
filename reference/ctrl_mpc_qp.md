# Inspect the MPC quadratic program

Return the QP data for given parameters, in OSQP form \\\min \frac12
z^\top H z + q^\top z\\, \\l \le A z \le u\\, together with index
helpers. This exposes the exact relation between control variables,
predicted states and solver matrices for verification.

## Usage

``` r
ctrl_mpc_qp(
  mpc,
  x,
  reference = NULL,
  u_prev = NULL,
  disturbance = NULL,
  k = 1L,
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

- output_offset:

  Optional known additive output offset (length `p`), for example an
  estimated output disturbance: output bounds then apply to \\C x + D
  u + \text{offset}\\.

## Value

A list with `H`, `q`, `A`, `l`, `u`, `constant` (objective constant) and
the index functions `xi(k)` and `ui(k)`.

## Examples

``` r
plant <- ctrl_ss(1, 1, 1, dt = 1)
qp <- ctrl_mpc_qp(ctrl_mpc(plant, 3, 1, 1), x = 1)
dim(qp$A)
#> [1] 4 7
```
