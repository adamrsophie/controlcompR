# Stability, controllability and matrix equations

``` r

library(controlcompR)
```

## Poles and stability

A continuous time model is asymptotically stable when every pole has a
negative real part; a discrete time model when every pole lies strictly
inside the unit circle. Poles on the boundary need more care: a simple
pole on the imaginary axis gives a marginally stable model, but a
repeated boundary pole with a nontrivial Jordan block grows without
bound.
[`ctrl_stability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_stability.md)
distinguishes these cases with a tolerance scaled by the size of `A`:

``` r

ctrl_stability(ctrl_ss(diag(c(-1, -2)), diag(2), diag(2)))$classification
#> [1] "asymptotically stable"
osc <- ctrl_ss(matrix(c(0, 1, -1, 0), 2, byrow = TRUE), c(0, 1), c(1, 0))
ctrl_stability(osc)$classification
#> [1] "marginally stable"
jordan <- ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0))
ctrl_stability(jordan)
#> Stability: unstable (continuous time)
#> Poles: 0+0i, 0+0i
#> Boundary poles (within tol = 4.44e-13): 0+0i, 0+0i
#> Stability margin: 0
```

[`ctrl_damp()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_damp.md)
lists natural frequencies and damping ratios (for discrete models
through the equivalent continuous poles \log(z)/T):

``` r

G <- ctrl_tf(9, c(1, 1.2, 9))
ctrl_damp(G)
#>             pole damping frequency time_constant
#> 1 -0.6+2.939388i     0.2         3      1.666667
#> 2 -0.6-2.939388i     0.2         3      1.666667
```

## Controllability and observability

The rank of the controllability matrix \[B, AB, \ldots, A^{n-1}B\] is a
poor numerical test: its columns grow like powers of `A`. The package
uses an orthogonal staircase reduction instead, which also isolates the
uncontrollable modes, and reports the Popov-Belevitch-Hautus test for
each eigenvalue so the *distance* to uncontrollability can be judged:

``` r

A <- diag(c(-1, 2))
res <- ctrl_controllability(A, c(1, 0))
res$rank
#> [1] 1
res$uncontrollable_modes
#> [1] 2+0i
res$stabilizable
#> [1] FALSE
res$pbh
#>   eigenvalue pbh_sigma_min controllable
#> 1       2+0i             0        FALSE
#> 2      -1+0i             1         TRUE
```

A model is *stabilizable* when every uncontrollable mode is stable and
*detectable* when every unobservable mode is stable. These weaker
properties are what LQR, Kalman filtering and MPC need:

``` r

ctrl_is_stabilizable(ctrl_ss(A, c(0, 1), diag(2)))
#> [1] TRUE
ctrl_observability(A, C = c(0, 1))$detectable
#> [1] TRUE
```

Rank decisions use tolerances. A nearly uncontrollable model is reported
as controllable with a tiny PBH singular value; the numbers matter more
than the yes or no answer:

``` r

near <- ctrl_controllability(diag(c(-1, -2)), c(1, 1e-9))
near$controllable
#> [1] TRUE
min(near$pbh$pbh_sigma_min)
#> [1] 7.071068e-10
```

## Lyapunov and Sylvester equations

[`ctrl_lyap()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lyap.md)
solves A X + X A^\top + Q = 0 (or A X A^\top - X + Q = 0) with the
Bartels-Stewart method on a real Schur form;
[`ctrl_sylvester()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_sylvester.md)
solves A X + X B = C. The solution of the Lyapunov equation with Q =
BB^\top is the controllability gramian:

``` r

sys <- ctrl_ss(matrix(c(-1, 2, 0, -3), 2, byrow = TRUE), c(1, 1), c(1, 0))
Wc <- ctrl_gram(sys, "c")
Wc
#>           [,1]      [,2]
#> [1,] 1.1666667 0.3333333
#> [2,] 0.3333333 0.1666667
#> attr(,"residual")
#> [1] 1.110223e-16
#> attr(,"separation")
#> [1] 0.6666667
max(abs(sys$A %*% Wc + Wc %*% t(sys$A) + sys$B %*% t(sys$B)))
#> [1] 2.220446e-16
```

The solvers check that the equation has a unique solution (no eigenvalue
pairs \lambda_i + \lambda_j = 0, or \lambda_i \lambda_j = 1 in discrete
time) and report the residual.

## Algebraic Riccati equations

[`ctrl_care()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_care.md)
and
[`ctrl_dare()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_care.md)
compute the *stabilizing* solutions of A^\top P + P A - (PB + S) R^{-1}
(B^\top P + S^\top) + Q = 0 and its discrete counterpart. The solvers
use the extended matrix pencil of size 2n + m, compress the input block
with a QR factorization (so `R` is never inverted) and compute an
ordered generalized Schur decomposition with LAPACK `dgges`. The
discrete solver does not need `A` to be invertible.

``` r

ric <- ctrl_care(sys$A, sys$B, diag(2), 1)
ric
#> Continuous algebraic Riccati equation: stabilizing solution
#> P =
#>        [,1]   [,2]
#> [1,] 0.3707 0.1379
#> [2,] 0.1379 0.2354
#> relative residual: 3.25e-16, symmetry error before symmetrization: 8.17e-17
#> closed-loop poles: -2.236+0i, -2.646+0i
```

The result carries the residual, the symmetry error before
symmetrization and the closed-loop poles, so the solution can be checked
rather than trusted. When no stabilizing solution exists, for example
because a mode on the imaginary axis is not reachable, the solver stops
with an explanation instead of returning a wrong matrix:

``` r

ctrl_care(osc$A, c(0, 0), diag(2), 1)
#> Error:
#> ! The computed Riccati solution is not stabilizing (closed-loop spectral abscissa 0); the problem is numerically too close to the boundary.
```
