# Conversions between model representations

- `ctrl_zpk2tf()` and `ctrl_tf2zpk()` convert between zero-pole-gain and
  polynomial forms (SISO).

- `ctrl_tf2ss()` realizes a proper transfer function in state space.
  SISO models use the controllable canonical form. MIMO models are
  realized element by element and then reduced with
  [`ctrl_minreal()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_minreal.md)
  when `minimal = TRUE`.

- `ctrl_ss2tf()` computes each SISO channel's transfer function. With
  `minimal = TRUE` (default) each channel is first reduced to a minimal
  realization, so uncontrollable or unobservable modes do not appear as
  canceling pole-zero pairs. The numerator of channel \\(i,j)\\ is
  computed from the identity \\c (sI - A)^{-1} b = \[\det(sI - A + bc) -
  \det(sI - A)\] / \det(sI - A)\\.

## Usage

``` r
ctrl_tf2ss(sys, minimal = TRUE, tol = NULL)

ctrl_ss2tf(sys, minimal = TRUE, tol = NULL)

ctrl_zpk2tf(sys)

ctrl_tf2zpk(sys)
```

## Arguments

- sys:

  A model.

- minimal:

  Reduce to a minimal realization (see details).

- tol:

  Tolerance passed to
  [`ctrl_minreal()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_minreal.md).

## Value

The converted model.

## Details

Root finding and characteristic polynomials are ill-conditioned for
high-order models; prefer state-space computations for large systems.

## Examples

``` r
G <- ctrl_tf(c(1, 2), c(1, 3, 2))
ss <- ctrl_tf2ss(G)
ss
#> State-space model (continuous time)
#> 1 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1
#> x1 -1
#> 
#> B =
#>        u1
#> x1 0.4472
#> 
#> C =
#>       x1
#> y1 2.236
#> 
#> D =
#>    u1
#> y1  0
ctrl_ss2tf(ss)
#> Transfer function (continuous time)
#> 
#>         1
#>   ---------
#>   s +     1
ctrl_tf2zpk(G)
#> Zero-pole-gain model (continuous time)
#>   gain: 1
#>   zeros: -2+0i
#>   poles: -2+0i, -1+0i
```
