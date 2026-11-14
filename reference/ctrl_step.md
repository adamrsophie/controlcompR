# Step, impulse and initial-condition responses

- `ctrl_step()`: response to a unit step (times `amplitude`) applied at
  time 0 to each selected input in turn, from zero initial state.

- `ctrl_impulse()`: impulse response. In continuous time this is the
  response to a Dirac impulse, computed as the free response from
  \\x(0^+) = B e_j\\. The Dirac component \\D\delta(t)\\ of the output
  cannot be represented on a time grid and is omitted (a warning is
  given when `D` is nonzero). In discrete time the input is the unit
  pulse (Kronecker delta) \\u_0 = 1\\, \\u_k = 0\\ for \\k \> 0\\; no
  `1/dt` scaling is applied, so a continuous impulse is *not* reproduced
  by a discrete pulse.

- `ctrl_initial()`: free response from the initial state `x0`.

## Usage

``` r
ctrl_step(sys, times = NULL, inputs = NULL, amplitude = 1)

ctrl_impulse(sys, times = NULL, inputs = NULL)

ctrl_initial(sys, x0, times = NULL)
```

## Arguments

- sys:

  A model.

- times:

  Optional time grid (sample instants for discrete models).

- inputs:

  Input channels (indices or names); default all.

- amplitude:

  Step amplitude.

- x0:

  Initial state for `ctrl_initial()`.

## Value

A `ctrl_response` with one case per input channel.

## Details

Default time grids are chosen from the pole locations.

## Examples

``` r
G <- ctrl_tf(4, c(1, 1.2, 4))
st <- ctrl_step(G)
ctrl_step_info(st)
#>   output    case rise_time settling_time overshoot undershoot     peak
#> 1     y1 step u1  0.661306      6.661667  37.10144          0 1.372324
#>   peak_time final_value settled
#> 1     1.645    1.000955    TRUE
plot(st)

ctrl_impulse(G, times = seq(0, 10, by = 0.05))
#> Time response (impulse, continuous time)
#>   201 time points from 0 to 10
#>   outputs: y1
```
