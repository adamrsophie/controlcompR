# LQG controller

Combine an LQR state feedback with a steady-state Kalman filter. Under
the separation principle (linear model, Gaussian noise with the given
covariances, and stabilizing LQR and Kalman designs) the closed-loop
poles are the union of the regulator poles, eig(A - BK), and the
estimator poles, eig(A - LC). LQG controllers have no guaranteed
robustness margins (Doyle, 1978); check margins of the actual loop.

## Usage

``` r
ctrl_lqg(
  sys,
  Q,
  R,
  Qn,
  Rn,
  N = NULL,
  Nn = NULL,
  G = NULL,
  form = c("current", "delayed")
)
```

## Arguments

- sys:

  Plant model.

- Q, R, N:

  LQR weights (see
  [`ctrl_lqr()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lqr.md)).

- Qn, Rn, Nn, G:

  Noise model (see
  [`ctrl_kalman()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman.md)).

- form:

  `"current"` or `"delayed"` (discrete models only).

## Value

An object of class `c("ctrl_lqg", "ctrl_controller")` with `controller`
(a `ctrl_ss` model from `y` to `u`), `lqr`, `kalman` and `form`.

## Details

The controller maps measurements `y` to inputs `u` (it contains the
minus sign of the feedback):

- continuous time and discrete `form = "delayed"`: \\\dot{\hat x} = (A -
  BK - LC + LDK)\hat x + L y\\, \\u = -K\hat x\\;

- discrete `form = "current"`: the input uses the measurement-updated
  estimate \\\hat x\_{k\|k}\\, \\u_k = -K \hat x\_{k\|k}\\, which reacts
  to \\y_k\\ without a one-sample delay.

## References

Doyle, J. C. (1978). Guaranteed margins for LQG regulators. *IEEE
Transactions on Automatic Control*, 23(4), 756-757.

## Examples

``` r
plant <- ctrl_c2d(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1),
                          c(1, 0)), 0.1)
lqg <- ctrl_lqg(plant, Q = diag(c(1, 0.1)), R = 0.1, Qn = diag(c(1e-4, 1e-3)),
                Rn = 1e-2)
cl <- ctrl_closed_loop(plant, lqg)
ctrl_poles(cl)
#> [1] 0.8677171+0.1006218i 0.8677171-0.1006218i 0.8675449+0.1003633i
#> [4] 0.8675449-0.1003633i
```
