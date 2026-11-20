# Stability analysis

Classify the internal stability of a linear time-invariant model from
its poles, with an explicit numerical tolerance.

## Usage

``` r
ctrl_stability(sys, tol = NULL, cluster_tol = NULL)
```

## Arguments

- sys:

  A model.

- tol:

  Boundary tolerance. Defaults to a multiple of machine precision times
  the norm of the balanced `A` (the state coordinates are first scaled
  to balance \\\[A\\ B; C\\ 0\]\\, so the result does not depend on
  state units, and a uniformly slow continuous model keeps its relative
  accuracy; in discrete time the scale is at least one). For a model
  whose time scales differ by about 1e13 or more, the slowest stable
  poles can fall inside this band and are classified as marginal; pass a
  smaller `tol` when such poles are known to be accurate.

- cluster_tol:

  Tolerance for grouping eigenvalues into clusters.

## Value

An object of class `ctrl_stability` with `classification`
(`"asymptotically stable"`, `"marginally stable"` or `"unstable"`),
`stable` (asymptotic stability), `poles`, `boundary_poles`,
`unstable_poles`, `stability_margin` (distance of the slowest pole to
the boundary; negative when unstable), the uniform band `tol`, the
per-pole bands `pole_tol` actually used and `cluster_tol`.

## Details

A continuous time model is asymptotically stable when every pole has
real part below `-tol`; a discrete time model when every pole has
modulus below `1 - tol`. Poles within `tol` of the stability boundary
(the imaginary axis or the unit circle) are *boundary poles*. If there
are boundary poles but no poles outside the boundary, the model is
*marginally stable* when every boundary eigenvalue is semisimple
(geometric multiplicity equals algebraic multiplicity, so trajectories
remain bounded) and *unstable* otherwise (a Jordan block on the boundary
produces polynomial growth). Multiplicities are numerical estimates:
eigenvalues are clustered within `cluster_tol` and the geometric
multiplicity uses a scale-aware SVD rank.

For transfer functions only the denominator roots are available, so the
Jordan structure test cannot be done; repeated boundary poles are then
reported as unstable (true for the input-output response).

## Examples

``` r
ctrl_stability(ctrl_ss(-1, 1, 1))
#> Stability: asymptotically stable (continuous time)
#> Poles: -1+0i
#> Stability margin: 1
# A double integrator is not marginally stable (Jordan block at 0)
ctrl_stability(ctrl_ss(matrix(c(0, 1, 0, 0), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#> Stability: unstable (continuous time)
#> Poles: 0+0i, 0+0i
#> Boundary poles (within tol = 4.44e-13): 0+0i, 0+0i
#> Stability margin: 0
# A harmonic oscillator is marginally stable
ctrl_stability(ctrl_ss(matrix(c(0, 1, -1, 0), 2, byrow = TRUE), c(0, 1), c(1, 0)))
#> Stability: marginally stable (continuous time)
#> Poles: 0+1i, 0-1i
#> Boundary poles (within tol = 4.44e-13): 0+1i, 0-1i
#> Stability margin: 0
```
