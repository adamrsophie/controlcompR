# Offset-free MPC estimator

Build a disturbance-augmented Kalman filter for offset-free tracking.
The plant model is augmented with a constant disturbance \\d\_{k+1} =
d_k\\ entering either the outputs (\\y = C x + D u + d\\,
`type = "output"`) or the inputs (\\x\_{k+1} = A x + B (u + d)\\,
`type = "input"`). With an observable augmented model and as many
disturbances as measured outputs, the closed loop has zero steady-state
tracking error for constant references and disturbances whenever it
converges, even under model mismatch.

## Usage

``` r
ctrl_mpc_estimator(mpc, type = c("output", "input"), Qn, Qd, Rn)
```

## Arguments

- mpc:

  A `ctrl_mpc` object.

- type:

  `"output"` or `"input"` disturbance model.

- Qn:

  Process noise covariance for the plant states.

- Qd:

  Covariance of the disturbance random walk.

- Rn:

  Measurement noise covariance.

## Value

An object of class `ctrl_mpc_estimator`.

## Details

At each sample the estimate \\(\hat x, \hat d)\\ is used to compute a
steady-state target \\(x_s, u_s)\\ from \$\$\begin{pmatrix} A - I & B \\
C & D \end{pmatrix} \begin{pmatrix} x_s \\ u_s \end{pmatrix} =
\begin{pmatrix} -B_d \hat d \\ r - C_d \hat d \end{pmatrix}\$\$ (least
squares when not square) and the MPC tracks it.

## Examples

``` r
plant <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.1)
mpc <- ctrl_mpc(plant, 10, Q = 0, Qy = 1, R = 0.01, u_min = -5, u_max = 5)
est <- ctrl_mpc_estimator(mpc, "output", Qn = 1e-4, Qd = 1e-3, Rn = 1e-4)
```
