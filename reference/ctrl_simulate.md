# Simulate a linear model

Compute the state and output trajectories of a linear time-invariant
model driven by an input signal.

## Usage

``` r
ctrl_simulate(
  sys,
  input = NULL,
  times = NULL,
  x0 = NULL,
  hold = c("zoh", "foh")
)
```

## Arguments

- sys:

  A model (transfer functions are realized with
  [`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)).

- input:

  Input signal (see details).

- times:

  Strictly increasing time points. For discrete models the default is
  `0:steps * dt` with `steps = nrow(input) - 1`.

- x0:

  Initial state (default zero).

- hold:

  `"zoh"` or `"foh"` input interpretation for continuous models.

## Value

A `ctrl_response` object with `time`, and arrays `x` (time x states x
1), `y` (time x outputs x 1) and `u` (time x inputs x 1), plus the
signal names. Use
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) for a
long data frame and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for a plot.

## Details

**Input format.** `input` is either `NULL` (zero input), a numeric
matrix with one row per time point and one column per input (a vector is
accepted for single-input models), or a function `function(t)` returning
the input vector at time `t` (it is evaluated on the time grid).

**Continuous time.** The simulation is exact for the chosen input
interpretation. With `hold = "zoh"` the input is held at `u(t_k)` on
\\\[t_k, t\_{k+1})\\; with `hold = "foh"` it is linearly interpolated
between `u(t_k)` and `u(t_{k+1})`. Each step uses the matrix exponential
of the corresponding block matrix (see
[`ctrl_c2d()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_c2d.md)),
so there is no integration error. The grid may be irregular; transition
matrices are cached for repeated step sizes.

**Discrete time.** The recurrence \\x\_{k+1} = A x_k + B u_k\\ is
evaluated exactly. `times` must be the sample instants \\t_0 + k\\dt\\;
`hold` is ignored.

Outputs follow \\y_k = C x_k + D u_k\\ at every time point.

## See also

[`ctrl_step()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md),
[`ctrl_impulse()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md),
[`ctrl_initial()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md)

## Examples

``` r
sys <- ctrl_ss(-1, 1, 1)
times <- seq(0, 5, by = 0.05)
r <- ctrl_simulate(sys, input = sin(times), times = times)
head(as.data.frame(r))
#>   time       case   kind signal       value
#> 1 0.00 simulation output     y1 0.000000000
#> 2 0.05 simulation output     y1 0.000000000
#> 3 0.10 simulation output     y1 0.002437513
#> 4 0.15 simulation output     y1 0.007187567
#> 5 0.20 simulation output     y1 0.014125209
#> 6 0.25 simulation output     y1 0.023125532
# Free response of x' = -x from x(0) = 1 equals exp(-t)
free <- ctrl_simulate(sys, times = times, x0 = 1)
max(abs(free$x[, 1, 1] - exp(-times)))
#> [1] 9.992007e-16
```
