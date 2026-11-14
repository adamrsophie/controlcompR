# Generalized plants and linear fractional transformations

A generalized plant has inputs `(w, u)` (exogenous inputs and controls)
and outputs `(z, y)` (performance outputs and measurements): \$\$\dot x
= A x + B_1 w + B_2 u, \quad z = C_1 x + D\_{11} w + D\_{12} u, \quad y
= C_2 x + D\_{21} w + D\_{22} u.\$\$ The last `ncon` inputs are controls
and the last `nmeas` outputs are measurements. `ctrl_lft()` closes the
loop `u = K y` and returns the map from `w` to `z`; it requires \\I -
D\_{22} D_K\\ to be invertible (well-posedness).

## Usage

``` r
ctrl_lft(P, K, nmeas, ncon)

ctrl_plant_partition(P, nmeas, ncon)
```

## Arguments

- P:

  Generalized plant (`ctrl_ss`).

- K:

  Controller from `y` to `u`.

- nmeas:

  Number of measurements `y`.

- ncon:

  Number of controls `u`.

## Value

`ctrl_lft()` returns the closed loop as a `ctrl_ss` model.
`ctrl_plant_partition()` returns the blocks as a list.

## Examples

``` r
P <- ctrl_ss(-1, cbind(1, 1), rbind(1, 1), matrix(0, 2, 2))
ctrl_lft(P, ctrl_gain(-2), nmeas = 1, ncon = 1)
#> State-space model (continuous time)
#> 1 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1
#> x1 -3
#> 
#> B =
#>    u1
#> x1  1
#> 
#> C =
#>    x1
#> y1  1
#> 
#> D =
#>    u1
#> y1  0
```
