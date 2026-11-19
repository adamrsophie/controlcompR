# Balanced realization and balanced truncation

`ctrl_balreal()` returns a balanced realization of an asymptotically
stable model, in which both Gramians equal \\\mathrm{diag}(\sigma_1,
\ldots, \sigma_r)\\ (the Hankel singular values); states with zero
Hankel singular value (to tolerance) are removed. `ctrl_balred()`
truncates the balanced realization to `order` states using the
square-root method of Tombs and Laub, which avoids forming the balancing
transformation from ill-conditioned Gramians.

## Usage

``` r
ctrl_balreal(sys, tol = 1e-12)

ctrl_balred(sys, order, method = c("truncate", "matchdc"))
```

## Arguments

- sys:

  A model.

- tol:

  Relative tolerance below which Hankel singular values are treated as
  zero in `ctrl_balreal()`.

- order:

  Number of states to keep.

- method:

  `"truncate"` or `"matchdc"`.

## Value

`ctrl_balreal()` returns a `ctrl_ss` model with attribute `hsv`.
`ctrl_balred()` returns the reduced `ctrl_ss` model with attributes
`hsv` (Hankel singular values of the stable part) and `error_bound`.

## Details

For an asymptotically stable model, balanced truncation preserves
stability and satisfies the a priori bound \$\$\\G - G_r\\\_\infty \le 2
\sum\_{i \> r} \sigma_i\$\$ (continuous and discrete time). With
`method = "matchdc"` the discarded states are residualized (singular
perturbation approximation), which matches the DC gain and satisfies the
same bound. The bound is checked against the error sampled on a
frequency grid, with the full model evaluated on a diagonally balanced
realization and its rounding error estimated from the conditioning of
\\zI - A\\. A sampled error above the bound by more than that accuracy
raises an error; otherwise `error_bound` is the larger of the a priori
bound and the sampled error plus its accuracy (relevant only when the
bound is at rounding level or attained, as when a single Hankel singular
value is discarded).

Unstable models are first split with
[`ctrl_stabsep()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_stabsep.md);
only the stable part is reduced and the unstable part is kept exactly,
so `order` must be at least the number of unstable poles.

## Examples

``` r
G <- ctrl_tf(c(1, 11, 30), c(1, 12, 47, 72, 36))
Gr <- ctrl_balred(G, 2)
attr(Gr, "error_bound")
#> [1] 0.005803089
ctrl_norm(ctrl_tf2ss(G) - Gr, "inf")
#> [1] 0.005803089
#> attr(,"frequency")
#> [1] 0
#> attr(,"bounds")
#> [1] 0.005803089 0.005803089
```
