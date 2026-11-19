# System norms

- `"2"`: the H2 norm. For an asymptotically stable continuous model it
  is \\\sqrt{\mathrm{tr}(C W_c C^\top)}\\ when `D = 0` and infinite
  otherwise; for a discrete model it is \\\sqrt{\mathrm{tr}(C W_c
  C^\top + D D^\top)}\\. The Gramian is computed after diagonal
  balancing; a realization whose Gramian is ill-conditioned for other
  reasons (nearly canceling input or output directions) loses relative
  accuracy, and
  [`ctrl_minreal()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_minreal.md)
  first helps.

- `"inf"`: the H-infinity norm \\\sup\_\omega
  \sigma\_{\max}(G(j\omega))\\ of an asymptotically stable model,
  computed with the Bruinsma-Steinbuch two-step algorithm. The lower
  bound comes from evaluations at candidate frequencies; each iteration
  tests a level \\\gamma\\ by checking the Hamiltonian matrix for
  imaginary-axis eigenvalues, which locate the frequency intervals where
  \\\sigma\_{\max} \> \gamma\\. The iteration stops when the relative
  gap between the bounds is below `tol`. Discrete models are mapped to
  continuous time with the norm-preserving bilinear transformation \\s =
  (z - 1)/(z + 1)\\. An upper bound is accepted only when no sampled
  gain on a verification grid (with local refinement) exceeds it;
  otherwise the search continues, and `bounds[2]` is `Inf` when no upper
  bound could be established.

## Usage

``` r
ctrl_norm(sys, type = c("2", "inf"), tol = 1e-10)
```

## Arguments

- sys:

  A model.

- type:

  `"2"` or `"inf"`.

- tol:

  Relative accuracy for the H-infinity norm.

## Value

A number; for `"inf"` with attribute `frequency` (the peak frequency in
rad per time unit) and `bounds`.

## Details

Unstable models have infinite H2 and H-infinity norms; `Inf` is
returned.

## References

Bruinsma, N. A. and Steinbuch, M. (1990). A fast algorithm to compute
the H-infinity norm of a transfer function matrix. *Systems and Control
Letters*, 14(4), 287-293.

## Examples

``` r
G <- ctrl_tf(1, c(1, 0.2, 1))
ctrl_norm(G, "inf")   # 1 / (2 * 0.1 * sqrt(1 - 0.01))
#> [1] 5.025189
#> attr(,"frequency")
#> [1] 0.9899495
#> attr(,"bounds")
#> [1] 5.025189 5.025189
ctrl_norm(ctrl_tf(1, c(1, 1)), "2")  # sqrt(1/2)
#> [1] 0.7071068
```
