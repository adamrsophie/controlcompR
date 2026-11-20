# Time-varying discrete Kalman filter and smoother

Run the Kalman filter on measured data for \\x\_{k+1} = A x_k + B u_k +
G w_k\\, \\y_k = C x_k + D u_k + v_k\\.

## Usage

``` r
ctrl_kalman_filter(
  sys,
  y,
  u = NULL,
  Qn,
  Rn,
  x0 = NULL,
  P0 = NULL,
  G = NULL,
  covariances = c("all", "variances", "final"),
  lag = NULL
)

ctrl_kalman_smoother(filter)
```

## Arguments

- sys:

  A discrete time `ctrl_ss` model.

- y:

  Measurements, a matrix with one row per sample (`T x p`).

- u:

  Inputs (`T x m`), or `NULL` for zero input.

- Qn, Rn:

  Process and measurement noise covariances.

- x0, P0:

  Prior mean and covariance at the first sample (defaults: zero and the
  identity).

- G:

  Process noise input matrix (default identity).

- covariances:

  `"all"` (default) stores the covariance, innovation covariance and
  gain of every sample; `"variances"` stores only the filtered
  variances, so memory grows with \\T n\\ instead of \\T n^2\\;
  `"final"` stores no per-sample covariance information, only the final
  covariances. The smoother needs `"all"`.

- lag:

  Optional fixed lag \\L \ge 0\\. When given, `lagged` holds the
  fixed-lag smoothed estimates: row `k` is the mean of \\x_k\\ given the
  measurements up to sample \\\min(k + L, T)\\. It needs only the last
  \\L + 1\\ smoother gains, so it works with every `covariances` mode.

  For online use, filter the data in chunks: pass `next_prediction` and
  `next_covariance` of one run as `x0` and `P0` of the next. The chunked
  estimates equal those of a single run over all the data.

- filter:

  A `ctrl_kalman_filter` result.

## Value

`ctrl_kalman_filter()` returns an object of class `ctrl_kalman_filter`
with `filtered` (`T x n`), `predicted` (`T x n`, row `k` is \\\hat
x\_{k\|k-1}\\), `variances` (`T x n`, the diagonal of \\P\_{k\|k}\\),
`P_filtered` and `P_predicted` (`n x n x T` arrays, `NULL` unless
`covariances = "all"`, as are `innovation_cov` and `gain`; `variances`
is `NULL` with `"final"`), `P_final` (the filtered covariance at the
last sample), `lagged` (`NULL` without `lag`), `innovations`,
`innovation_cov` (`p x p x T`), `gain` (`n x p x T`), `nis` (normalized
innovation squared), `loglik` (Gaussian log-likelihood of the
measurements) and `next_prediction`, `next_covariance` (the prediction
for the sample after the data). `ctrl_kalman_smoother()` returns a list
with `smoothed` and `P_smoothed`.

## Details

**Conventions.** `x0` and `P0` are the prior mean and covariance of the
state at the first sample, \\\hat x\_{1\|0}\\ and \\P\_{1\|0}\\. Each
step first *updates* with \\y_k\\ and then *predicts* to \\k + 1\\:
\$\$e_k = y_k - C \hat x\_{k\|k-1} - D u_k, \quad S_k = C P\_{k\|k-1}
C^\top + R_n, \quad K_k = P\_{k\|k-1} C^\top S_k^{-1},\$\$ \$\$\hat
x\_{k\|k} = \hat x\_{k\|k-1} + K_k e_k, \quad \hat x\_{k+1\|k} = A \hat
x\_{k\|k} + B u_k.\$\$ The covariance update uses the Joseph form
\\P\_{k\|k} = (I - K_k C) P\_{k\|k-1} (I - K_k C)^\top + K_k R_n
K_k^\top\\, which keeps it symmetric positive semidefinite in floating
point. Rows of `y` that are entirely `NA` are treated as missing
measurements (the update is skipped); partially missing rows, infinite
values and non-numeric data are errors.

`ctrl_kalman_smoother()` applies the Rauch-Tung-Striebel backward pass
to the filter output.

## Examples

``` r
set.seed(1)
sys <- ctrl_ss(matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE), c(0.005, 0.1),
               c(1, 0), dt = 0.1)
x <- c(0, 1)
y <- numeric(50)
for (k in 1:50) {
  y[k] <- x[1] + rnorm(1, sd = 0.1)
  x <- drop(sys$A %*% x) + rnorm(2, sd = c(0.001, 0.01))
}
kf <- ctrl_kalman_filter(sys, y, Qn = diag(c(1e-6, 1e-4)), Rn = 0.01)
tail(kf$filtered, 3)
#>             x1        x2
#> [48,] 4.680705 0.9807787
#> [49,] 4.764152 0.9704786
#> [50,] 4.863031 0.9717680
sm <- ctrl_kalman_smoother(kf)
```
