# PID controllers

Define a PID controller with the features needed in practice: derivative
filtering, setpoint weighting, output saturation, anti-windup and
bumpless transfer. The control law (parallel form) is \$\$u = k_p
(b\\r - y) + k_i \int (r - y)\\dt + k_d \frac{d}{dt} (c\\r - y)\_f,\$\$
where the derivative acts on the filtered signal \\(c r - y)\_f = (c r -
y) / (T_f s + 1)\\ and `b`, `c` are setpoint weights (`c = 0` avoids
derivative kicks on setpoint steps). With `form = "standard"` the gains
are given as \\k_p\\, \\T_i\\ and \\T_d\\ with \\k_i = k_p / T_i\\ and
\\k_d = k_p T_d\\.

## Usage

``` r
ctrl_pid(
  kp = 1,
  ki = 0,
  kd = 0,
  Tf = NULL,
  derivative_filter = NULL,
  b = 1,
  c = 0,
  u_min = -Inf,
  u_max = Inf,
  antiwindup = c("backcalc", "clamp", "none"),
  kt = NULL,
  form = c("parallel", "standard"),
  dt = NULL,
  integral_method = c("forward", "backward", "tustin")
)
```

## Arguments

- kp:

  Proportional gain.

- ki:

  Integral gain (parallel form) or integral time `Ti` (standard form;
  `Inf` disables integral action).

- kd:

  Derivative gain (parallel form) or derivative time `Td` (standard
  form).

- Tf:

  Derivative filter time constant.

- derivative_filter:

  Alternative to `Tf`: the ratio `N` with `Tf = Td / N`.

- b, c:

  Setpoint weights for the proportional and derivative terms.

- u_min, u_max:

  Output saturation limits.

- antiwindup:

  `"backcalc"`, `"clamp"` or `"none"`.

- kt:

  Back-calculation tracking gain.

- form:

  `"parallel"` or `"standard"`.

- dt:

  Sample interval of the discrete implementation (`NULL` for a
  continuous description; required for simulation).

- integral_method:

  Discretization of the integral term.

## Value

An object of class `c("ctrl_pid", "ctrl_controller")`.

## Details

The derivative filter time constant `Tf` is required whenever \\k_d \ne
0\\: an unfiltered derivative is improper and amplifies measurement
noise without bound. A common choice is \\T_f = T_d / N\\ with \\N\\
between 2 and 20 (`derivative_filter = N`).

**Discrete implementation** (used by
[`ctrl_simulate_pid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_simulate_pid.md)
and
[`ctrl_pid_update()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_update.md))
with sample interval `dt`: backward Euler for the derivative filter,
forward Euler for the integral by default
(`integral_method = "forward"`; `"backward"` and `"tustin"` are also
available), and saturation to `[u_min, u_max]`. Anti-windup options:
`"backcalc"` (back-calculation, the integrator is driven by \\k_t
(u\_{sat} - u)\\ with tracking gain `kt`, default \\1/\sqrt{T_i T_d}\\
or \\1/T_i\\), `"clamp"` (conditional integration: the integrator is
frozen while saturated in the direction of the error) or `"none"`.
Bumpless transfer from manual to automatic mode initializes the
integrator so that the first automatic output equals the last manual
output.

## See also

[`ctrl_pid_tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_tf.md),
[`ctrl_simulate_pid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_simulate_pid.md),
[`ctrl_pid_tune()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_tune.md)

## Examples

``` r
pid <- ctrl_pid(kp = 2, ki = 1, kd = 0.5, derivative_filter = 10,
                u_min = -1, u_max = 1, dt = 0.01)
pid
#> PID controller (discrete, dt = 0.01)
#>   kp = 2, ki = 1, kd = 0.5, Tf = 0.025
#>   setpoint weights b = 1, c = 0
#>   saturation [-1, 1], anti-windup: backcalc
ctrl_pid_tf(pid)
#> Transfer function (discrete time, dt = 0.01)
#> 
#>   16.29 z^2 - 31.99 z + 15.71
#>   ---------------------------
#>     z^2 - 1.714 z + 0.7143
```
