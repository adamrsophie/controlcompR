# Sensitivity functions

For a loop transfer function `L` closed with unity negative feedback,
\\S = (I + L)^{-1}\\ and \\T = L (I + L)^{-1} = I - S\\. With a plant
`P` and controller `C` (loop \\L = P C\\), `ctrl_gangof4()` also returns
\\S P = (I + PC)^{-1} P\\ (input disturbance to output) and \\C S\\
(output noise to control signal).

## Usage

``` r
ctrl_sensitivity(L)

ctrl_gangof4(P, C)
```

## Arguments

- L:

  Loop transfer function.

- P, C:

  Plant and controller.

## Value

A list of models `S` and `T` (and `PS`, `CS`).

## Examples

``` r
L <- ctrl_tf(10, c(1, 3, 3, 1))
sens <- ctrl_sensitivity(L)
ctrl_norm(sens$S, "inf")  # peak sensitivity Ms
#> [1] Inf
```
