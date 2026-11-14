# Stochastic closed-loop simulation of an LQG controller

Simulate a discrete time plant with Gaussian process and measurement
noise under an LQG controller, optionally over several independent
realizations, and report the true and estimated state trajectories.
Summaries over realizations (mean squared estimation error, mean
quadratic cost) avoid judging a stochastic design from one favorable
trajectory.

## Usage

``` r
ctrl_simulate_lqg(
  plant,
  lqg,
  steps,
  x0,
  Qn = NULL,
  Rn = NULL,
  replications = 1L,
  seed = NULL
)
```

## Arguments

- plant:

  Discrete time plant (`ctrl_ss`); it may differ from the design model
  to study model mismatch, but must have the same dimensions.

- lqg:

  A `ctrl_lqg` object.

- steps:

  Number of samples.

- x0:

  Initial true state.

- Qn, Rn:

  Covariances used to *generate* the noise (default: those of the
  design).

- replications:

  Number of independent realizations.

- seed:

  Optional seed (the previous RNG state is restored).

## Value

A list of class `ctrl_lqg_simulation` with arrays `x`, `xhat`, `u`, `y`
(samples x signals x replications) and a `summary` data frame with the
per-replication RMS estimation error and average stage cost.

## Examples

``` r
plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                          c(1, 0)), 0.1)
lqg <- ctrl_lqg(plant, Q = diag(c(1, 0.1)), R = 0.1, Qn = diag(c(1e-4, 1e-3)),
                Rn = 1e-2)
sim <- ctrl_simulate_lqg(plant, lqg, steps = 100, x0 = c(1, 0),
                         replications = 5, seed = 42)
sim$summary
#>   replication rms_estimation_error average_cost
#> 1           1            0.1861161    0.1422334
#> 2           2            0.1964052    0.1568127
#> 3           3            0.2053522    0.1576309
#> 4           4            0.2048822    0.1837231
#> 5           5            0.1718826    0.1342959
```
