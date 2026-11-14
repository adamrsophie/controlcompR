# Natural frequencies and damping ratios

For a continuous time pole \\\lambda\\, the natural frequency is
\\\omega_n = \|\lambda\|\\ and the damping ratio is \\\zeta =
-\mathrm{Re}(\lambda)/\|\lambda\|\\. Discrete time poles are first
mapped to the equivalent continuous time pole \\\log(\lambda)/dt\\;
poles at the origin of the z plane have no continuous equivalent and get
`NA` values.

## Usage

``` r
ctrl_damp(sys)
```

## Arguments

- sys:

  A model.

## Value

A data frame with `pole`, `damping`, `frequency` (rad per time unit) and
`time_constant` (`1 / (zeta * wn)`, `Inf` when undamped).

## Examples

``` r
ctrl_damp(ctrl_tf(1, c(1, 0.4, 4)))
#>             pole damping frequency time_constant
#> 1 -0.2+1.989975i     0.1         2             5
#> 2 -0.2-1.989975i     0.1         2             5
```
