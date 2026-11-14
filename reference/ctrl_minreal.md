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

  A model. The default tolerance first tries
  \\\sqrt{\epsilon}\max(\\A\\, \\B\\, \\C\\)\\ (as in MATLAB), which
  removes modes that cancel to roughly half the working precision, as
  needed for realizations built from polynomial coefficients. The
  reduced model is then compared with the original at frequencies spread
  over the range of the poles; when they disagree (a badly scaled model
  whose small couplings were mistaken for zeros), the reduction is
  repeated with the much tighter staircase tolerance \\10 n^2 \epsilon
  \\\[A, B\]\\\\, and an error is raised if even that result does not
  match.

- tol:

  Absolute rank tolerance for the staircase reductions. When given, it
  is used as is and the frequency check is skipped.

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
