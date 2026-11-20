# Mixed-sensitivity generalized plant

Build the generalized plant for the mixed-sensitivity problem
\$\$\left\\ \begin{pmatrix} W_1 S \\ W_2 K S \\ W_3 T \end{pmatrix}
\right\\\$\$ with \\S = (I + G K)^{-1}\\ and \\T = G K S\\. The plant
has input `(w, u)` with `w` the reference (or output disturbance) and
output `(z1, z2, z3, e)` with \\e = w - G u\\; the controller uses `e`
and the loop convention is \\u = K e\\. Any weight may be `NULL`
(omitted).

## Usage

``` r
ctrl_augw(G, W1 = NULL, W2 = NULL, W3 = NULL)
```

## Arguments

- G:

  Plant.

- W1, W2, W3:

  Weights on S, KS and T (`ctrl_tf`, `ctrl_ss` or numbers).

## Value

A `ctrl_ss` generalized plant with `nmeas = ncon =` the plant
dimensions.

## Details

For use with
[`ctrl_hinfsyn()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_hinfsyn.md)
(which requires \\D\_{11} = 0\\), `W1` must be strictly proper and `W2`
must have full rank feedthrough (so \\D\_{12}\\ has full column rank).
`W3 G` must be proper. Plant poles on the imaginary axis (integrators)
are not driven by `w` and make the H-infinity Riccati equations
singular; shift them slightly into the left half plane for the design
and verify the controller on the true plant.

## Examples

``` r
G <- ctrl_tf(1, c(1, 3, 2))
W1 <- ctrl_tf(1, c(1, 0.01))
P <- ctrl_augw(G, W1 = W1, W2 = 0.1)
syn <- ctrl_hinfsyn(P, nmeas = 1, ncon = 1)
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
syn$gamma
#> [1] 0.6896671
```
