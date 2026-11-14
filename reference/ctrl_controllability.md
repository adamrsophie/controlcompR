# Controllability and stabilizability analysis

Determine the controllable subspace with an orthogonal staircase
reduction, which is far better conditioned than the rank of the
controllability matrix. The uncontrollable modes are the eigenvalues of
the uncontrollable block. The model is *stabilizable* when every
uncontrollable mode is asymptotically stable; full controllability is
not required.

## Usage

``` r
ctrl_controllability(sys, B = NULL, tol = NULL, discrete = NULL)
```

## Arguments

- sys:

  A `ctrl_ss` model (or the matrix `A`).

- B:

  Input matrix when `sys` is a matrix.

- tol:

  Rank tolerance.

- discrete:

  Time domain when `sys` is a matrix (`FALSE` by default).

## Value

A list of class `ctrl_controllability` with `rank`, `n`, `controllable`,
`stabilizable`, `uncontrollable_modes`, `basis` (an orthonormal basis of
the controllable subspace), `transform` (the orthogonal staircase
transformation), `pbh` (table of PBH tests) and `ctrb_singular_values`.

## Details

Rank decisions use the tolerance `tol`, scaled by the norms of `A` and
`B`. A rank computed under a tolerance is a numerical statement: a
nearly uncontrollable system can be reported as controllable with a tiny
PBH singular value. Inspect `pbh$pbh_sigma_min` (the smallest singular
value of \\\[A - \lambda I, B\]\\ for each eigenvalue) to judge the
distance to uncontrollability.

## Examples

``` r
A <- diag(c(-1, 2))
B <- c(1, 0)
res <- ctrl_controllability(A, B)
res$controllable
#> [1] FALSE
res$stabilizable  # the unstable mode at 2 is uncontrollable
#> [1] FALSE
```
