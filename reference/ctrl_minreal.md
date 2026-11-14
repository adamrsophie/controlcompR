# Minimal realization

Remove uncontrollable and unobservable states with two orthogonal
staircase reductions (see
[`ctrl_controllability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_controllability.md)).
The input-output behavior is preserved up to the rank tolerance.

## Usage

``` r
ctrl_minreal(sys, tol = NULL)
```

## Arguments

- sys:

  A model.

- tol:

  Absolute rank tolerance for the staircase reductions. The default,
  \\\sqrt{\epsilon}\max(\\A\\, \\B\\, \\C\\)\\ as in MATLAB, removes
  modes that cancel to roughly half the working precision, which is
  needed for realizations built from polynomial coefficients.

## Value

A `ctrl_ss` model with a minimal number of states.

## Examples

``` r
# Two copies of the same first-order system: one state is redundant
sys <- ctrl_ss(diag(c(-1, -1)), c(1, 1), c(1, 1))
ctrl_dims(ctrl_minreal(sys))
#>  states  inputs outputs 
#>       1       1       1 
```
