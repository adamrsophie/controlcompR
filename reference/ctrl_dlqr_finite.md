# Finite-horizon discrete LQR

Solve the finite-horizon problem \$\$\min \sum\_{k=0}^{N-1} (x_k^\top Q
x_k + u_k^\top R u_k) + x_N^\top P_N x_N\$\$ with the backward Riccati
recursion \$\$P_k = Q + A^\top P\_{k+1} A - A^\top P\_{k+1} B (R +
B^\top P\_{k+1} B)^{-1} B^\top P\_{k+1} A.\$\$ The optimal input is the
time-varying feedback \\u_k = -K_k x_k\\. This is the exact reference
for unconstrained linear MPC with the same weights and horizon.

## Usage

``` r
ctrl_dlqr_finite(sys, Q, R, horizon, P_terminal = NULL)
```

## Arguments

- sys:

  A discrete time model.

- Q, R:

  Weights (scalars and vectors are expanded to diagonals).

- horizon:

  Number of steps `N`.

- P_terminal:

  Terminal weight \\P_N\\ (defaults to `Q`).

## Value

A list with `K` (list of `N` gains, `K[[1]]` applies at `k = 0`) and `P`
(list of `N + 1` cost matrices, `P[[1]]` at `k = 0`).

## Examples

``` r
sys <- ctrl_ss(1, 1, 1, dt = 1)
fin <- ctrl_dlqr_finite(sys, 1, 1, horizon = 30)
fin$K[[1]]  # approaches the infinite-horizon gain
#>          [,1]
#> [1,] 0.618034
ctrl_lqr(sys, 1, 1)$K
#>          x1
#> u1 0.618034
```
