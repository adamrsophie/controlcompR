# Robust performance under multiplicative output uncertainty

Test whether the performance specification \\\\W_p S_p\\\_\infty \< 1\\
holds for every perturbed plant \\G_p = (I + W_u \Delta) G\\ with
\\\\\Delta\\\_\infty \le 1\\, where \\S_p = (I + G_p K)^{-1}\\. Since
\\S_p = S (I + W_u \Delta T)^{-1}\\, the condition \$\$\sup\_\omega
\left(\|W_p(j\omega)\|\\ \bar\sigma(S(j\omega)) + \|W_u(j\omega)\|\\
\bar\sigma(T(j\omega))\right) \< 1\$\$ together with nominal internal
stability is sufficient. For SISO loops it is also necessary (the
classical robust performance condition). For MIMO loops with the scalar
weights used here it can be conservative; it is not a structured
singular value computation.

## Usage

``` r
ctrl_robust_performance(G, K, Wp, Wu, n = 2000L)
```

## Arguments

- G:

  Nominal plant.

- K:

  Controller (negative feedback).

- Wp:

  Scalar performance weight (SISO system or number).

- Wu:

  Scalar uncertainty weight (SISO system or number).

- n:

  Number of grid frequencies.

## Value

A list with `robust_performance` (sampled peak below one), `certified`
(`TRUE` when the rigorous upper bound is also below one), `upper_bound`,
`nominally_stable`, `robust_stability` (\\\sup \|W_u\| \bar\sigma(T) \<
1\\), `nominal_performance` (\\\sup \|W_p\| \bar\sigma(S) \< 1\\),
`peak`, `frequency` (where the peak occurs; `Inf` for the infinite
frequency endpoint), `peaks` (the sampled sum, performance and
robustness peaks), `criterion_exact` (`TRUE` for SISO loops, where the
condition is also necessary) and `curve` (a data frame with `frequency`
and `value`).

## Details

Each supremum (performance, robustness and their sum) is sampled on a
dense frequency grid that includes the endpoints (zero and infinite
frequency for continuous models, zero and the Nyquist frequency for
discrete models), and the largest local maxima are refined with a
one-dimensional search. Sampled peaks are lower bounds on the true
suprema, so `robust_performance` is decided on samples. A pass is
`certified` when the rigorous upper bound \\\\W_p S\\\_\infty + \\W_u
T\\\_\infty\\ from
[`ctrl_norm()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_norm.md)
is below one as well; it can be conservative, so an uncertified pass is
not a failure.

## Examples

``` r
G <- ctrl_tf(1, c(1, 1))
K <- ctrl_tf(c(2, 1), c(1, 0))
Wp <- ctrl_tf(c(0.5, 0.5), c(1, 0.01))
Wu <- ctrl_tf(c(1, 0.2), c(0.1, 1))
rp <- ctrl_robust_performance(G, K, Wp, Wu)
rp$peak
#> [1] 2.051606
```
