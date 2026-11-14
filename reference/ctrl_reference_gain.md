# Reference (feedforward) gain for state feedback

For the law \\u = -K x + N_r r\\, compute \\N_r\\ so that the closed
loop has unit DC gain from `r` to the outputs `y` (square case, least
squares otherwise). This gives zero steady-state error for constant
references only for the exact nominal model; use integral action or
offset-free MPC for robustness to model error.

## Usage

``` r
ctrl_reference_gain(plant, K)
```

## Arguments

- plant:

  Plant model.

- K:

  Gain matrix or `ctrl_lqr` object.

## Value

The `m x p` matrix \\N_r\\.

## Examples

``` r
plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
lqr <- ctrl_lqr(plant, diag(2), 1)
ctrl_reference_gain(plant, lqr)
#>          [,1]
#> v_u1 2.236068
```
