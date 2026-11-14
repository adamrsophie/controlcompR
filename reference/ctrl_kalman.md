# Steady-state Kalman filter

For the model \$\$x\_{k+1} = A x_k + B u_k + G w_k, \quad y_k = C x_k +
D u_k + v_k\$\$ (or its continuous time analogue) with \\E\[w w^\top\] =
Q_n\\, \\E\[v v^\top\] = R_n\\, \\E\[w v^\top\] = N_n\\, compute the
steady-state Kalman gains from the filtering Riccati equation (the dual
of the LQR equation).

## Usage

``` r
ctrl_kalman(sys, Qn, Rn, Nn = NULL, G = NULL)
```

## Arguments

- sys:

  A `ctrl_ss` model.

- Qn:

  Process noise covariance (size `ncol(G)`).

- Rn:

  Measurement noise covariance (`p x p`).

- Nn:

  Optional cross covariance (`ncol(G) x p`).

- G:

  Process noise input matrix; defaults to the identity (noise on every
  state).

## Value

An object of class `c("ctrl_kalman", "ctrl_estimator")` with `L`, `M`
(discrete only), `P`, `Z` (discrete only), `poles` (eigenvalues of \\A -
L C\\), and `system` (the predictor as a `ctrl_ss` model with inputs `u`
then `y`).

## Details

Discrete time conventions: `P` is the steady-state *prior* covariance
\\P\_{k\|k-1}\\; `M` is the measurement-update (innovation) gain \\M = P
C^\top (C P C^\top + R_n)^{-1}\\ used in \\\hat x\_{k\|k} = \hat
x\_{k\|k-1} + M (y_k - C\hat x\_{k\|k-1} - D u_k)\\; `L` is the one-step
predictor gain \\L = (A P C^\top + G N_n)(C P C^\top + R_n)^{-1}\\ used
in \\\hat x\_{k+1\|k} = A \hat x\_{k\|k-1} + B u_k + L(y_k - C \hat
x\_{k\|k-1} - D u_k)\\; `Z` is the posterior covariance \\P\_{k\|k}\\
(for \\N_n = 0\\). In continuous time, `L` is the Kalman-Bucy gain \\(P
C^\top + G N_n) R_n^{-1}\\.

`Rn` must be positive definite: perfect (noise-free) measurements make
the filtering problem singular and are not silently regularized. `Qn`
may be singular (including zero) as long as the Riccati equation has a
stabilizing solution; this requires \\(A, C)\\ detectable and no
unreachable noise modes on the stability boundary.

## See also

[`ctrl_kalman_filter()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman_filter.md)
for the time-varying filter.

## Examples

``` r
sys <- ctrl_ss(matrix(c(1, 0.1, 0, 1), 2, byrow = TRUE), c(0.005, 0.1),
               c(1, 0), dt = 0.1)
kf <- ctrl_kalman(sys, Qn = diag(c(1e-4, 1e-3)), Rn = 0.01)
kf$M
#>           [,1]
#> [1,] 0.2372931
#> [2,] 0.2761715
```
