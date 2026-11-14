# Root locus

Closed-loop poles of \\1 + k L(s) = 0\\ for a SISO loop as the gain `k`
varies, computed as the roots of \\d(s) + k\\n(s)\\. Branches are
tracked by matching each root to the nearest root at the previous gain.

## Usage

``` r
ctrl_rlocus(L, gains = NULL)
```

## Arguments

- L:

  SISO loop transfer function.

- gains:

  Optional nonnegative gains; by default a logarithmic grid chosen so
  that the poles approach the open-loop zeros or asymptotes.

## Value

An object of class `ctrl_rlocus` with `gains` and a complex matrix
`roots` (gains x branches).

## Examples

``` r
rl <- ctrl_rlocus(ctrl_tf(1, c(1, 3, 2, 0)))
dim(rl$roots)
#> [1] 401   3
```
