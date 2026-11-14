# Series, parallel and feedback connections

- `ctrl_series(sys1, sys2)`: the output of `sys1` drives `sys2`, so the
  overall transfer function is \\G_2 G_1\\ (note the order). The
  operator form is `sys2 * sys1`.

- `ctrl_parallel(sys1, sys2)`: \\G_1 + G_2\\ (same inputs, summed
  outputs); operator form `sys1 + sys2`.

- `ctrl_feedback(sys1, sys2, sign = -1)`: closes the loop \\u_1 = r +
  \mathrm{sign} \cdot y_2\\, \\u_2 = y_1\\, with output \\y_1\\. The
  default `sign = -1` is negative feedback, giving \\(I + G_1 G_2)^{-1}
  G_1\\.

- `ctrl_append(...)`: block-diagonal (parallel, unconnected) stacking.

## Usage

``` r
ctrl_series(sys1, sys2)

ctrl_parallel(sys1, sys2)

ctrl_feedback(sys1, sys2 = 1, sign = -1)

ctrl_append(...)
```

## Arguments

- sys1, sys2:

  Models or numeric gains.

- sign:

  Feedback sign, `-1` (negative) or `+1` (positive).

- ...:

  Models to append.

## Value

A model.

## Details

When both models have direct feedthrough the feedback loop contains an
algebraic loop; it is resolved by solving with \\I - \mathrm{sign}\\ D_1
D_2\\, and an error is raised when this matrix is singular or
ill-conditioned (the interconnection is then not well posed).

Models must share the time domain; a static gain adopts the time domain
of the other model. SISO transfer functions are combined with polynomial
arithmetic (without cancellation); all other combinations are carried
out in state space.

## Examples

``` r
G <- ctrl_tf(1, c(1, 1))
C <- ctrl_tf(c(2, 1), c(1, 0))
L <- ctrl_series(C, G)        # G * C
T <- ctrl_feedback(L, 1)      # unity negative feedback
ctrl_poles(T)
#> [1] -0.381966+0i -2.618034+0i
```
