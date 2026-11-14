# Discretize a continuous time model

Convert a continuous time model to discrete time with sample interval
`dt`.

## Usage

``` r
ctrl_c2d(
  sys,
  dt,
  method = c("zoh", "foh", "tustin", "matched"),
  prewarp = NULL
)
```

## Arguments

- sys:

  A continuous time model.

- dt:

  Positive sample interval.

- method:

  One of `"zoh"`, `"foh"`, `"tustin"` or `"matched"`.

- prewarp:

  Optional prewarping frequency for `"tustin"`.

## Value

A discrete time model of the same class as `sys` (`"matched"` returns a
transfer function).

## Details

- `"zoh"` (default): exact zero-order hold. The input is held constant
  between samples and \\A_d = e^{A T}\\, \\B_d = \int_0^T e^{A\tau}
  B\\d\tau\\, both read from the matrix exponential of the block matrix
  \\\begin{pmatrix} A & B \\ 0 & 0 \end{pmatrix} T\\. This is valid for
  singular `A`; the formula \\A^{-1}(A_d - I)B\\ is never used.

- `"foh"`: exact first-order (triangle) hold, for inputs that are
  linearly interpolated between samples. The discrete state is \\\xi_k =
  x_k - \Gamma_2 u_k\\ (not the physical state), so a physical initial
  state \\x_0\\ corresponds to \\\xi_0 = x_0 - \Gamma_2 u_0\\; the
  matrix \\\Gamma_2\\ is stored in the `foh_state_map` attribute of the
  result.

- `"tustin"`: bilinear transformation \\s =
  \frac{2}{T}\frac{z-1}{z+1}\\, optionally prewarped so that the
  frequency response matches exactly at `prewarp` rad per time unit.

- `"matched"`: matched pole-zero mapping \\z = e^{sT}\\ for SISO models,
  with zeros at infinity mapped to \\z = -1\\ (all but one, so the
  result stays strictly proper when the original is) and the gain
  matched at DC (or at `s = j/T` when there is a pole or zero at 0).

Matrix exponentials use
[`expm::expm()`](https://rdrr.io/pkg/expm/man/expm.html) (scaling and
squaring with Pade approximation).

## Examples

``` r
sys <- ctrl_ss(-1, 1, 1)
sysd <- ctrl_c2d(sys, 0.1)
sysd$A  # exp(-0.1)
#>           [,1]
#> [1,] 0.9048374
ctrl_c2d(ctrl_tf(1, c(1, 1, 0)), 0.5, method = "tustin")
#> Transfer function (discrete time, dt = 0.5)
#> 
#>    0.05 z^2 +   0.1 z +  0.05
#>   ---------------------------
#>      z^2 -   1.6 z +   0.6
```
