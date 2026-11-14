# Observability and detectability analysis

The dual of
[`ctrl_controllability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_controllability.md):
the analysis is applied to \\(A^\top, C^\top)\\. The model is
*detectable* when every unobservable mode is asymptotically stable.

## Usage

``` r
ctrl_observability(sys, C = NULL, tol = NULL, discrete = NULL)
```

## Arguments

- sys:

  A `ctrl_ss` model (or the matrix `A`).

- C:

  Output matrix when `sys` is a matrix.

- tol:

  Rank tolerance.

- discrete:

  Time domain when `sys` is a matrix (`FALSE` by default).

## Value

A list of class `ctrl_controllability` with `rank`, `observable`,
`detectable`, `unobservable_modes`, `basis` (orthonormal basis of the
observable subspace, that is the orthogonal complement of the
unobservable subspace), `pbh` and `ctrb_singular_values` (singular
values of the observability matrix).

## Examples

``` r
ctrl_observability(diag(c(-1, -2)), C = c(1, 0))$detectable
#> [1] TRUE
```
