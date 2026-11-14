# Gain, phase and delay margins

Classical stability margins of a SISO loop transfer function `L` (the
loop is closed with unity negative feedback). Gain crossover frequencies
(\\\|L\| = 1\\) and phase crossover frequencies (\\L\\ real and
negative) are computed as roots of polynomials in \\\omega\\ (or on the
unit circle for discrete models), then polished with a bracketing root
finder, so they do not depend on a frequency grid.

## Usage

``` r
ctrl_margin(L)
```

## Arguments

- L:

  A SISO loop transfer function (`ctrl_tf`, `ctrl_ss` or `ctrl_zpk`).

## Value

A list of class `ctrl_margin` with `gain_margin` (absolute),
`gain_margin_db`, `phase_margin` (degrees), `delay_margin`,
`phase_crossover` and `gain_crossover` frequencies, the tables `all_gm`
and `all_pm`, and `closed_loop_stable`.

## Details

The gain margin at a phase crossover \\\omega\_{pc}\\ is
\\1/\|L(j\omega\_{pc})\|\\; the phase margin at a gain crossover
\\\omega\_{gc}\\ is \\180 + \angle L(j\omega\_{gc})\\ degrees, wrapped
to \\(-180, 180\]\\; the delay margin is the smallest \\(\mathrm{PM}\_i
\bmod 360)/\omega\_{gc,i}\\ over all gain crossovers (radians over rad
per time unit), the smallest added delay that brings a crossover onto
the critical point. The reported gain and phase margins are the most
critical ones (gain margin closest to 1 on a log scale, smallest phase
margin). `closed_loop_stable` is computed on the given realization, so
unstable hidden modes of a state-space loop are not cancelled away.
Margins describe distance to instability only if the nominal closed loop
is stable, which is reported in `closed_loop_stable`. Margins are SISO
quantities; for MIMO loops use loop-at-a-time disk margins.

## Examples

``` r
L <- ctrl_tf(10, c(1, 3, 3, 1))
ctrl_margin(L)
#> Stability margins (unity negative feedback)
#>   gain margin: 0.8 (-1.938 dB) at w = 1.732
#>   phase margin: -7.033 deg at w = 1.908
#>   delay margin: 3.228
#>   nominal closed loop stable: FALSE
```
