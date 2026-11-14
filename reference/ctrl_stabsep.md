# Stable and unstable decomposition

Split a model into `G = G_stable + G_unstable` with an ordered real
Schur form followed by a Sylvester-equation block diagonalization. The
feed through `D` is assigned to the stable part. Poles on the stability
boundary (within `tol`) make the split ill-defined and cause an error.

## Usage

``` r
ctrl_stabsep(sys, tol = NULL)
```

## Arguments

- sys:

  A model.

- tol:

  Boundary tolerance.

## Value

A list with `stable` and `unstable` `ctrl_ss` models.

## Examples

``` r
sys <- ctrl_ss(diag(c(-1, 2, -3)), c(1, 1, 1), c(1, 1, 1))
parts <- ctrl_stabsep(sys)
ctrl_poles(parts$unstable)
#> [1] 2+0i
```
