# Closed-loop MPC simulation

Run the receding-horizon loop: at each sample, measure (or estimate) the
state, solve the MPC problem, apply the first input to the plant and
advance. The simulation plant may differ from the prediction model
(model mismatch) and may include input disturbances and measurement
noise.

## Usage

``` r
ctrl_simulate_mpc(
  mpc,
  x0,
  steps,
  reference = NULL,
  plant = NULL,
  disturbance = NULL,
  noise = NULL,
  estimator = NULL,
  u0 = NULL,
  on_failure = c("stop", "hold", "shifted"),
  warm_start = TRUE
)
```

## Arguments

- mpc:

  A `ctrl_mpc` object.

- x0:

  Initial plant state.

- steps:

  Number of samples.

- reference:

  List with `x`, `u` and/or `y` references (see
  [`ctrl_mpc_qp()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mpc_qp.md));
  with an estimator, `y` is the output setpoint.

- plant:

  Simulation plant (default: the prediction model).

- disturbance:

  Input disturbance added to the plant input: a vector with one value
  per input (constant), a matrix with one row per sample, or a function
  of the sample index. For a single input, a vector with one value per
  sample is also accepted. The plant sees the effective input \\u_k +
  d_k\\ in both \\x\_{k+1} = A x_k + B(u_k + d_k)\\ and \\y_k = C x_k +
  D(u_k + d_k)\\. The measurement used to compute \\u_k\\ is taken
  before \\u_k\\ is applied, so it sees the input and disturbance held
  over the previous sample (zero before the first).

- noise:

  Measurement noise on the outputs (same formats).

- estimator:

  Optional `ctrl_mpc_estimator` for output feedback with offset-free
  tracking; without it the true state is measured.

- u0:

  Input applied before the first sample (for rate constraints).

- on_failure:

  `"stop"`, `"hold"` or `"shifted"`.

- warm_start:

  Warm start each solve with the shifted previous plan.

## Value

An object of class `ctrl_mpc_simulation` with `time`, `x` (true states,
`(steps+1) x n`), `y`, `u`, `xhat`, `dhat`, a `log` data frame (status,
iterations, solve time, objective, max violation, fallback) and summary
counts.

## Details

**Failure handling.** If a QP solve fails, `on_failure` decides what
happens; the choice is always recorded in the `fallback` column:
`"stop"` (default) raises an error with the solver status, `"hold"`
applies the previous input (clipped to the input bounds) and `"shifted"`
applies the next move of the last successful plan (then holds its last
move). No input is ever applied silently.

## Examples

``` r
plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                          c(1, 0)), 0.1)
mpc <- ctrl_mpc(plant, horizon = 20, Q = diag(c(10, 1)), R = 0.1,
                u_min = -1, u_max = 1, x_min = c(-Inf, -0.6),
                x_max = c(Inf, 0.6))
sim <- ctrl_simulate_mpc(mpc, x0 = c(2, 0), steps = 80)
summary(sim)
#> MPC simulation: 80 samples, 0 failed solves (0%)
#>   solver statuses: solved = 80
#>   solve time: mean 0.000478 s, max 0.00143 s; mean iterations 406.6
#>   max input bound violation: 9.43e-11
#>   max state bound violation (true plant): 0
#>   max slack: 0
```
