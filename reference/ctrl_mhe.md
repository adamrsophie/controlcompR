# Moving horizon estimation for linear models

Estimate the state of \\x\_{k+1} = A x_k + B u_k + G w_k\\, \\y_k = C
x_k + D u_k + v_k\\ by solving, at each sample \\T\\, the finite-horizon
problem over the window \\s = T - N + 1, \ldots, T\\: \$\$\min\_{x_s,
\ldots, x_T,\\ w_s, \ldots, w\_{T-1}} \\x_s - \bar x_s\\^2\_{P_s^{-1}} +
\sum\_{j=s}^{T-1} \\w_j\\^2\_{Q_n^{-1}} + \sum\_{j=s}^{T} \\y_j - C
x_j - D u_j\\^2\_{R_n^{-1}}\$\$ subject to the model equations and
optional bounds on the states and the process noise. The estimate is
\\\hat x\_{T\|T} = x_T\\.

## Usage

``` r
ctrl_mhe(
  sys,
  horizon,
  Qn,
  Rn,
  G = NULL,
  x0 = NULL,
  P0 = NULL,
  x_min = -Inf,
  x_max = Inf,
  w_min = -Inf,
  w_max = Inf,
  solver_settings = list(),
  accuracy = c("high", "standard", "fast")
)

ctrl_mhe_run(mhe, y, u = NULL, warm_start = TRUE)
```

## Arguments

- sys:

  Discrete time `ctrl_ss` model.

- horizon:

  Window length `N` (number of measurements in the window).

- Qn, Rn:

  Process and measurement noise covariances (positive definite).

- G:

  Process noise input matrix (default identity).

- x0, P0:

  Prior mean and covariance at the first sample.

- x_min, x_max:

  Bounds on the states.

- w_min, w_max:

  Bounds on the process noise.

- solver_settings:

  OSQP settings; they override the `accuracy` profile.

- accuracy:

  Solver accuracy profile, as in
  [`ctrl_mpc()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mpc.md)
  (`"high"` uses tolerances 1e-10).

- mhe:

  A `ctrl_mhe` object.

- y:

  Measurements (`T x p`, `NA` rows are missing samples).

- u:

  Inputs (`T x m`) or `NULL`.

- warm_start:

  Start each window's QP from the previous window's solution, shifted by
  one sample (the new last state is the model prediction and the new
  last disturbance is zero). This changes only the number of OSQP
  iterations, not the solution.

## Value

An object of class `ctrl_mhe`; run it with `ctrl_mhe_run()`.

A list of class `ctrl_mhe_result` with `estimates` (`T x n`, row `k` is
\\\hat x\_{k\|k}\\), `status` (OSQP status per sample), `objective`,
`iterations` and `solve_time`.

## Details

**Arrival cost.** The first term summarizes the data before the window.
The prior mean is the prediction from the previous MHE estimate, \\\bar
x_s = A \hat x\_{s-1\|s-1} + B u\_{s-1}\\, and \\P_s\\ follows the
Kalman filter covariance recursion. Without active constraints this
makes MHE coincide with the Kalman filter (which the tests verify); with
active constraints MHE differs from, and can improve on, the Kalman
filter, at the price of one QP per sample, so its cost grows with the
horizon. MHE is worth it when states or noises have hard physical bounds
(concentrations, saturations) that the Kalman filter would violate.

Each window QP is solved by OSQP. Once the window is full, the solver
workspace is reused from sample to sample while the sparsity pattern of
the QP does not change; only its numbers are updated.

## Examples

``` r
sys <- ctrl_ss(0.9, 1, 1, dt = 1)
mhe <- ctrl_mhe(sys, horizon = 5, Qn = 0.01, Rn = 0.1, x_min = 0)
est <- ctrl_mhe_run(mhe, y = c(0.2, -0.1, 0.4, 0.3, 0.1, 0.5))
est$estimates
#>              x1
#> [1,] 0.18181818
#> [2,] 0.04356436
#> [3,] 0.15438123
#> [4,] 0.18145126
#> [5,] 0.14818606
#> [6,] 0.21654918
```
