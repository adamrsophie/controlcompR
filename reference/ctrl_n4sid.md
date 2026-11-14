# Subspace state-space identification (PO-MOESP)

Estimate a discrete time state-space model \\x\_{k+1} = A x_k + B u_k\\,
\\y_k = C x_k + D u_k + v_k\\ from input-output data with the
past-output MOESP subspace method (Verhaegen, 1994). Future outputs are
projected on future inputs and instrumented by past inputs and outputs
(an LQ factorization of stacked block Hankel matrices); the singular
values of the projection reveal the order, the extended observability
matrix gives `A` and `C` (shift invariance), and `B`, `D` and the
initial state are then estimated by linear least squares on the
simulated output.

## Usage

``` r
ctrl_n4sid(data, order = NULL, horizon = NULL, feedthrough = FALSE)
```

## Arguments

- data:

  A `ctrl_iddata` object.

- order:

  Model order; `NULL` chooses the order at the largest gap of the log
  singular values.

- horizon:

  Number of block rows `i` (must exceed the order).

- feedthrough:

  Estimate `D` (default `FALSE`, so `D = 0`).

## Value

An object of class `c("ctrl_n4sid", "ctrl_idmodel")` with `sys` (a
`ctrl_ss`), `x0`, `singular_values`, `order` and `residuals`.

## Details

The instruments make the estimates consistent for white or colored
output noise in open loop. Closed-loop data (inputs correlated with past
noise) bias the estimates; this is not corrected. The model has no noise
model (output-error form).

## References

Verhaegen, M. (1994). Identification of the deterministic part of MIMO
state space models given in innovations form from input-output data.
*Automatica*, 30(1), 61-74.

## Examples

``` r
set.seed(2)
true <- ctrl_ss(matrix(c(0.7, 0.2, -0.2, 0.7), 2), c(1, 0), c(1, 1), dt = 1)
u <- rnorm(400)
y <- ctrl_simulate(true, input = u, times = 0:399)$y[, 1, 1] + rnorm(400, sd = 0.02)
m <- ctrl_n4sid(ctrl_iddata(y, u), order = 2)
ctrl_poles(m$sys)
#> [1] 0.6989143+0.1993426i 0.6989143-0.1993426i
```
