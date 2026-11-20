# Robust stability under unstructured uncertainty

Test robust stability of the feedback loop of plant `G` and controller
`K` (negative feedback, \\u = -K y\\) against a family of unstructured,
stable, norm-bounded perturbations \\\Delta\\ with \\\\\Delta\\\_\infty
\le 1\\:

## Usage

``` r
ctrl_robust_stability(
  G,
  K,
  W,
  type = c("multiplicative_output", "multiplicative_input", "additive")
)
```

## Arguments

- G:

  Nominal plant.

- K:

  Controller (negative feedback).

- W:

  Uncertainty weight (system or number).

- type:

  Uncertainty structure.

## Value

A list with `robustly_stable` (decided on `upper_bound`),
`nominally_stable`, `peak` (the H-infinity norm tested), `upper_bound`
(the certified upper bound on that norm from
[`ctrl_norm()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_norm.md)),
`margin` (`1 / peak`) and `type`.

## Details

- `"multiplicative_output"`: \\G\_\Delta = (I + W \Delta) G\\; robustly
  stable iff \\\\T W\\\_\infty \< 1\\, \\T = G K (I + G K)^{-1}\\;

- `"multiplicative_input"`: \\G\_\Delta = G (I + W \Delta)\\; iff
  \\\\T_i W\\\_\infty \< 1\\, \\T_i = K G (I + K G)^{-1}\\;

- `"additive"`: \\G\_\Delta = G + W \Delta\\; iff \\\\K S W\\\_\infty \<
  1\\, \\S = (I + G K)^{-1}\\.

**Guarantee.** By the small-gain theorem, if the nominal loop is
internally stable and the test value is below one, the loop is stable
for *every* perturbation in the family (an exact condition for
unstructured complex \\\Delta\\; the perturbed plants must keep the same
number of unstable poles as `G` for the additive and multiplicative
forms). The `margin` is the largest scaling of the uncertainty that is
tolerated. The result says nothing about structured or real parametric
uncertainty, for which it can be conservative, and nothing about robust
*performance* (see
[`ctrl_robust_performance()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_robust_performance.md)).

## Examples

``` r
G <- ctrl_tf(1, c(1, 1))
K <- ctrl_tf(10, 1)
W <- ctrl_tf(c(1, 0.2), c(0.1, 1))
ctrl_robust_stability(G, K, W)
#> $robustly_stable
#> [1] FALSE
#> 
#> $nominally_stable
#> [1] TRUE
#> 
#> $peak
#> [1] 4.762771
#> 
#> $upper_bound
#> [1] 4.762771
#> 
#> $margin
#> [1] 0.2099618
#> 
#> $type
#> [1] "multiplicative_output"
#> 
```
