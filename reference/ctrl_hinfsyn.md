# H-infinity output-feedback synthesis (experimental)

Compute a stabilizing controller with closed-loop \\\\T\_{zw}\\\_\infty
\< \gamma\\ for a continuous time generalized plant (see
[`ctrl_lft()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lft.md)),
and optionally search for the smallest achievable \\\gamma\\ by
bisection.

## Usage

``` r
ctrl_hinfsyn(P, nmeas, ncon, gamma = NULL, tol = 0.001, gamma_range = NULL)
```

## Arguments

- P:

  Generalized plant.

- nmeas, ncon:

  Numbers of measurements and controls.

- gamma:

  Target level; `NULL` searches for the optimum.

- tol:

  Relative bisection tolerance.

- gamma_range:

  Optional `c(lower, upper)` search interval.

## Value

An object of class `c("ctrl_synthesis", "ctrl_controller")` with `K`,
`closed_loop`, `gamma` (the level used), `norm` (certified closed-loop
H-infinity norm), `X`, `Y`, `history` (bisection log) and `method`.

## Details

**Scope.** Continuous time only, with \\D\_{11} = 0\\, \\D\_{12}\\ full
column rank, \\D\_{21}\\ full row rank, \\(A, B_2)\\ stabilizable and
\\(C_2, A)\\ detectable; arbitrary \\D\_{22}\\ is handled with a loop
shift. Plants with \\D\_{11} \ne 0\\ (for example mixed-sensitivity
designs with a biproper sensitivity weight) are rejected, as are
discrete time plants.

**Algorithm.** The plant is first normalized with orthogonal
transformations of `z` and `w` and invertible scalings of `u` and `y` so
that \\D\_{12} = \[0; I\]\\ and \\D\_{21} = \[0, I\]\\ (these leave the
H-infinity norm unchanged). For a given \\\gamma\\, the two H-infinity
Riccati equations of Glover and Doyle (1988) are solved with
[`ctrl_care()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_care.md)
on indefinite weights; \\\gamma\\ is feasible when both stabilizing
solutions exist, \\X\_\infty \succeq 0\\, \\Y\_\infty \succeq 0\\ and
\\\rho(X\_\infty Y\_\infty) \< \gamma^2\\. The central controller is
then \\A_K = A + \gamma^{-2} B_1 B_1^\top X\_\infty + B_2 F\_\infty +
Z\_\infty L\_\infty (C_2 + \gamma^{-2} D\_{21} B_1^\top X\_\infty)\\,
\\B_K = -Z\_\infty L\_\infty\\, \\C_K = F\_\infty\\, \\D_K = 0\\, with
\\F\_\infty = -(D\_{12}^\top C_1 + B_2^\top X\_\infty)\\, \\L\_\infty =
-(B_1 D\_{21}^\top + Y\_\infty C_2^\top)\\ and \\Z\_\infty = (I -
\gamma^{-2} Y\_\infty X\_\infty)^{-1}\\.

**Certification.** Every returned controller is checked independently:
the closed loop must be internally stable and its H-infinity norm,
computed with
[`ctrl_norm()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_norm.md)
(Bruinsma-Steinbuch, a different algorithm), must be below \\\gamma\\.
If the check fails, an error is raised instead of returning an
uncertified controller. The bisection stops at relative tolerance `tol`
and returns the controller for the smallest feasible \\\gamma\\ found,
which is slightly suboptimal by design: the central controller becomes
ill-conditioned at the optimum. When the H2 optimal controller already
gives a zero closed-loop norm, the optimum zero is attained and that
controller is returned with `gamma = 0`.

## References

Glover, K. and Doyle, J. C. (1988). State-space formulae for all
stabilizing controllers that satisfy an H-infinity norm bound and
relations to risk sensitivity. *Systems and Control Letters*, 11,
167-172.

## Examples

``` r
A <- matrix(c(0, 1, 0, 0), 2, byrow = TRUE)
P <- ctrl_ss(A, cbind(diag(2), 0, c(0, 1)),
             rbind(c(1, 0), c(0, 0), c(1, 0)),
             rbind(c(0, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0)))
hinf <- ctrl_hinfsyn(P, nmeas = 1, ncon = 1)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
c(gamma = hinf$gamma, certified_norm = hinf$norm)
#>          gamma certified_norm 
#>       2.896528       2.896523 
```
