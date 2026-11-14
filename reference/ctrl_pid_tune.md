# PID tuning rules

- `"ziegler_nichols"`: the classical frequency-response rule from the
  ultimate gain \\K_u\\ (the gain margin of the plant) and the ultimate
  period \\P_u = 2\pi/\omega\_{180}\\: \\k_p = 0.6 K_u\\, \\T_i =
  P_u/2\\, \\T_d = P_u/8\\ (PID) or \\k_p = 0.45 K_u\\, \\T_i = P_u /
  1.2\\ (PI). It gives aggressive, poorly damped loops and is included
  as a reference point.

- `"simc"`: Skogestad's SIMC rules for a first-order plus dead-time
  model \\K e^{-\theta s}/(\tau s + 1)\\ with closed-loop time constant
  \\\tau_c\\ (default \\\theta\\): \\k_p = \tau / (K(\tau_c +
  \theta))\\, \\T_i = \min(\tau, 4(\tau_c + \theta))\\ (PI).

- `"amigo"`: the AMIGO PI rule of Astrom and Hagglund for the same
  model.

## Usage

``` r
ctrl_pid_tune(
  method = c("ziegler_nichols", "simc", "amigo"),
  plant = NULL,
  fopdt = NULL,
  type = c("PID", "PI"),
  tau_c = NULL,
  ...
)
```

## Arguments

- method:

  Tuning rule.

- plant:

  Plant model (Ziegler-Nichols).

- fopdt:

  Numeric vector `c(K, tau, theta)`.

- type:

  `"PID"` or `"PI"` (Ziegler-Nichols only; SIMC and AMIGO give PI
  controllers).

- tau_c:

  SIMC closed-loop time constant.

- ...:

  Further arguments for
  [`ctrl_pid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid.md)
  (for example `dt`, `u_min`, `u_max`, `derivative_filter`).

## Value

A `ctrl_pid` object.

## Details

Supply `plant` for Ziegler-Nichols, or `fopdt = c(K, tau, theta)` (for
example from
[`ctrl_fit_fopdt()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_fit_fopdt.md))
for SIMC and AMIGO. The rules are heuristics; verify the tuned loop by
simulation and margins.

## Examples

``` r
plant <- ctrl_tf(1, c(1, 3, 3, 1))
ctrl_pid_tune("ziegler_nichols", plant = plant, derivative_filter = 10)
#> PID controller (continuous description)
#>   kp = 4.8, ki = 2.646, kd = 2.177, Tf = 0.04534
#>   setpoint weights b = 1, c = 0
#>   saturation [-Inf, Inf], anti-windup: backcalc
ctrl_pid_tune("simc", fopdt = c(1, 2, 0.5))
#> PID controller (continuous description)
#>   kp = 2, ki = 1, kd = 0, Tf = 0
#>   setpoint weights b = 1, c = 0
#>   saturation [-Inf, Inf], anti-windup: backcalc
```
