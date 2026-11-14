# Luenberger observer

Design an observer gain `L` by placing the eigenvalues of \\A - L C\\
with
[`ctrl_place()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_place.md)
applied to the dual pair \\(A^\top, C^\top)\\. The observer is
\$\$\dot{\hat x} = A \hat x + B u + L (y - C \hat x - D u)\$\$ (or the
corresponding difference equation in discrete time, a one-step
predictor).

## Usage

``` r
ctrl_observer(sys, poles)
```

## Arguments

- sys:

  A `ctrl_ss` model (observable).

- poles:

  Desired observer poles.

## Value

An object of class `c("ctrl_observer", "ctrl_estimator")` with `L`,
`poles` and `system` (the observer as a `ctrl_ss` model with inputs `u`
then `y` and the state estimates as outputs).

## Examples

``` r
plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
obs <- ctrl_observer(plant, c(-10, -11))
obs$L
#>      [,1]
#> [1,]   18
#> [2,]   54
#> attr(,"poles")
#> [1] -11+0i -10+0i
#> attr(,"cond_X")
#> [1] 114.0088
#> attr(,"max_error")
#> [1] 1.385558e-13
```
