# Linear quadratic regulator

Design the state feedback \\u = -K x\\ that minimizes, over inputs that
stabilize the closed loop, \$\$J = \int_0^\infty \left(x^\top Q x +
u^\top R u + 2 x^\top N u \right) dt\$\$ for a continuous time model, or
\$\$J = \sum\_{k=0}^\infty \left(x_k^\top Q x_k + u_k^\top R u_k + 2
x_k^\top N u_k\right)\$\$ for a discrete time model. The gain is \\K =
R^{-1}(B^\top P + N^\top)\\ (continuous) or \\K = (R + B^\top P
B)^{-1}(B^\top P A + N^\top)\\ (discrete), where `P` is the stabilizing
solution of
[`ctrl_care()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_care.md)
or
[`ctrl_dare()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_care.md).

## Usage

``` r
ctrl_lqr(sys, Q, R, N = NULL)
```

## Arguments

- sys:

  A model.

- Q:

  State weight (`n x n`). A scalar or vector is expanded to a diagonal
  matrix.

- R:

  Input weight (`m x m`). A scalar or vector is expanded to a diagonal
  matrix.

- N:

  Optional cross weight (`n x m`).

## Value

An object of class `c("ctrl_lqr", "ctrl_controller")` with `K`, `P`,
`poles` (closed-loop eigenvalues), `residual` (Riccati residual), the
weights and the plant.

## Details

**Sign convention.** The control law is \\u = -K x\\, so the closed loop
is \\A - B K\\.

Prerequisites are checked: `R` must be positive definite, the joint
weight \\\[Q, N; N^\top, R\]\\ positive semidefinite and \\(A, B)\\
stabilizable. Without stabilization as a requirement the minimum can be
smaller: an unstable mode that the cost does not see (for example with
`Q = 0`) is cheapest left alone. The two coincide when the cost is
detectable. The optimality statement applies to the linear,
unconstrained, exactly known model; LQR is not optimal for nonlinear,
constrained or uncertain plants, and the guaranteed margins of
continuous time LQR (infinite gain margin, 60 degree phase margin per
loop when \\N = 0\\ and `R` is diagonal) hold only with full state
measurement.

## See also

[`ctrl_closed_loop()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_closed_loop.md),
[`ctrl_dlqr_finite()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_dlqr_finite.md),
[`ctrl_lqg()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lqg.md)

## Examples

``` r
plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
controller <- ctrl_lqr(plant, Q = diag(c(10, 1)), R = 1)
controller$K
#>          x1        x2
#> u1 1.741657 0.6719633
closed <- ctrl_closed_loop(plant, controller)
ctrl_stability(closed)
#> Stability: asymptotically stable (continuous time)
#> Poles: -1.836+0.609i, -1.836-0.609i
#> Stability margin: 1.836
```
