# Frequency response

Evaluate \\G(j\omega)\\ for a continuous time model or
\\G(e^{j\omega\\dt})\\ for a discrete time model. For discrete models
frequencies above the Nyquist frequency \\\pi/dt\\ alias and produce a
warning.

## Usage

``` r
ctrl_freqresp(sys, w = NULL)
```

## Arguments

- sys:

  A model.

- w:

  Frequencies in rad per time unit. If `NULL`, a logarithmic grid is
  chosen from the pole and zero locations (see
  [`ctrl_freq_grid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_freq_grid.md)).

## Value

A list of class `ctrl_freqresp` with `w` and `response`, a complex array
of dimension outputs x inputs x frequencies.

## Details

State-space models are balanced (diagonal scaling) and evaluated by
solving \\(j\omega I - A) X = B\\ at each frequency (no matrix inverse).
Transfer functions are evaluated as ratios of their numerator and
denominator polynomials (Horner's rule), which can lose accuracy for
high order polynomials; convert to state space with
[`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
for those.

## Examples

``` r
G <- ctrl_tf(1, c(1, 1))
fr <- ctrl_freqresp(G, c(0.1, 1, 10))
Mod(fr$response[1, 1, ])
#> [1] 0.99503719 0.70710678 0.09950372
```
