# Closed-loop systems

Build the closed loop of a plant with a controller.

## Usage

``` r
ctrl_closed_loop(plant, controller, include_control = TRUE)
```

## Arguments

- plant:

  Plant model.

- controller:

  A `ctrl_lqr` object, a gain matrix `K`, a `ctrl_lqg` object, or a
  `ctrl_ss` / `ctrl_tf` error-feedback controller.

- include_control:

  Append the control signal to the outputs.

## Value

A `ctrl_ss` model.

## Details

- With an LQR or other static state feedback (`ctrl_lqr`, or a matrix
  `K`), the control law is \\u = -K x + v\\, where `v` is an external
  input. The closed loop is \\\dot x = (A - BK) x + B v\\, \\y = (C -
  DK) x + D v\\, with outputs `y` followed (when
  `include_control = TRUE`) by the control signal `u`.

- With an LQG controller (see
  [`ctrl_lqg()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lqg.md))
  the controller is a dynamic system from the measurement to the input,
  connected in positive feedback `u = K_c(y)` (the controller already
  contains the minus sign), with an additive input disturbance `v` at
  the plant input. States are the plant states followed by the
  controller states.

- With a dynamic `ctrl_ss` or `ctrl_tf` controller `Cs` acting on the
  error \\r - y\\, the loop is the unity negative feedback \\y = (I + G
  C_s)^{-1} G C_s r\\.

## Examples

``` r
plant <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
cl <- ctrl_closed_loop(plant, ctrl_lqr(plant, diag(2), 1))
ctrl_poles(cl)
#> [1] -1.000000+0i -2.236068+0i
```
