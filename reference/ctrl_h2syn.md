# H2 optimal output-feedback synthesis

Compute the controller minimizing the H2 norm of the closed loop
\\T\_{zw}\\ for a continuous time generalized plant (see
[`ctrl_lft()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lft.md)).

## Usage

``` r
ctrl_h2syn(P, nmeas, ncon)
```

## Arguments

- P:

  Generalized plant.

- nmeas, ncon:

  Numbers of measurements and controls.

## Value

An object of class `c("ctrl_synthesis", "ctrl_controller")` with `K`
(controller), `closed_loop`, `norm` (closed-loop H2 norm), `X`, `Y` and
`method`.

## Details

**Assumptions** (checked; violations raise errors): \\D\_{11} = 0\\ (so
the closed-loop H2 norm is finite), \\D\_{12}\\ full column rank,
\\D\_{21}\\ full row rank, \\(A, B_2)\\ stabilizable, \\(C_2, A)\\
detectable, and no invariant zeros of the \\(A, B_2, C_1, D\_{12})\\ and
\\(A, B_1, C_2, D\_{21})\\ subsystems on the imaginary axis (detected as
Riccati failures).

**Algorithm.** Two Riccati equations (Zhou, Doyle and Glover, 1996,
Theorem 14.7): \\X_2\\ = `ctrl_care(A, B2, C1'C1, D12'D12, C1'D12)`,
\\F_2 = -(D\_{12}^\top D\_{12})^{-1}(B_2^\top X_2 + D\_{12}^\top C_1)\\;
\\Y_2\\ = `ctrl_care(A', C2', B1 B1', D21 D21', B1 D21')`, \\L_2 = -(Y_2
C_2^\top + B_1 D\_{21}^\top)(D\_{21} D\_{21}^\top)^{-1}\\; the
controller is \\A_K = A + B_2 F_2 + L_2 C_2\\, \\B_K = -L_2\\, \\C_K =
F_2\\, \\D_K = 0\\ (with a loop shift when \\D\_{22} \ne 0\\). The
closed-loop H2 norm is computed independently with
[`ctrl_norm()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_norm.md)
and returned.

## References

Zhou, K., Doyle, J. C. and Glover, K. (1996). *Robust and Optimal
Control*. Prentice Hall.

## Examples

``` r
# LQG as an H2 problem: z = [x1; u], w = [process noise; sensor noise]
A <- matrix(c(0, 1, 0, 0), 2, byrow = TRUE)
P <- ctrl_ss(A, cbind(diag(2), 0, c(0, 1)),
             rbind(c(1, 0), c(0, 0), c(1, 0)),
             rbind(c(0, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0)))
h2 <- ctrl_h2syn(P, nmeas = 1, ncon = 1)
h2$norm
#> [1] 3.29439
```
