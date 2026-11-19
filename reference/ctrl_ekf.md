# Extended and unscented Kalman filters

Nonlinear state estimation for a discrete time `ctrl_nlsys` with
additive noise, \\x\_{k+1} = f(k, x_k, u_k) + w_k\\, \\y_k = h(k, x_k,
u_k) + v_k\\, using the conventions of
[`ctrl_kalman_filter()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman_filter.md)
(`x0`, `P0` are the prior at the first sample; each step updates with
\\y_k\\ and then predicts; rows of `y` that are entirely `NA` skip the
update). Row `k` of the data is passed to `f` and `h` as the sample
index `k - 1`, the 0-based convention of
[`ctrl_nl_simulate()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nl_simulate.md).

## Usage

``` r
ctrl_ekf(sys, y, u = NULL, Qn, Rn, x0 = NULL, P0 = NULL)

ctrl_ukf(
  sys,
  y,
  u = NULL,
  Qn,
  Rn,
  x0 = NULL,
  P0 = NULL,
  alpha = 1,
  beta = 2,
  kappa = 0
)
```

## Arguments

- sys:

  Discrete time `ctrl_nlsys`.

- y:

  Measurements (`T x p`).

- u:

  Inputs (`T x m`) or `NULL`.

- Qn, Rn:

  Process and measurement noise covariances.

- x0, P0:

  Prior mean and covariance.

- alpha, beta, kappa:

  Unscented transform parameters.

## Value

A list of class `c("ctrl_nl_filter", "ctrl_kalman_filter")` with
`filtered`, `P_filtered`, `innovations`, `nis` and `loglik`. It has no
predicted covariances, so
[`ctrl_kalman_smoother()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman_filter.md)
does not accept it.

## Details

- `ctrl_ekf()` linearizes `f` and `h` at the current estimate with
  central finite differences and uses the Joseph-form covariance update.

- `ctrl_ukf()` propagates \\2n + 1\\ sigma points (scaled unscented
  transform with parameters `alpha`, `beta`, `kappa`); it needs no
  Jacobians and is exact for linear models.

Both are approximations for nonlinear models: their covariances are not
guaranteed to be consistent, and they can diverge for strong
nonlinearities or poor initialization. Check the normalized innovation
squared (`nis`), whose mean should be close to the number of outputs.

## Examples

``` r
f <- function(k, x, u, p) c(x[1] + 0.1 * x[2], x[2] - 0.1 * sin(x[1]))
h <- function(k, x, u, p) x[1]
sys <- ctrl_nlsys(f, h, nstates = 2, ninputs = 0, dt = 0.1)
set.seed(1)
x <- c(1, 0)
y <- numeric(100)
for (k in 1:100) {
  y[k] <- x[1] + rnorm(1, sd = 0.05)
  x <- f(k, x, NULL, NULL) + rnorm(2, sd = 0.001)
}
ekf <- ctrl_ekf(sys, y, Qn = 1e-6, Rn = 0.0025, x0 = c(0.8, 0))
ukf <- ctrl_ukf(sys, y, Qn = 1e-6, Rn = 0.0025, x0 = c(0.8, 0))
tail(cbind(ekf$filtered, ukf$filtered), 2)
#>               x1         x2        x1         x2
#>  [99,] -1.276753 -0.7291844 -1.276617 -0.7291842
#> [100,] -1.348630 -0.6332541 -1.348497 -0.6332617
```
