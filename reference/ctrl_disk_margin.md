# Disk margins

The disk margin of a SISO loop `L` (unity negative feedback) is the
largest \\\alpha\\ such that the loop stays stable for every complex
gain perturbation \\f\\ in the disk \$\$D(\alpha, \sigma) =
\left\\\frac{1 + \frac{1-\sigma}{2}\delta} {1 -
\frac{1+\sigma}{2}\delta} : \|\delta\| \le \alpha\right\\,\$\$ where the
skew \\\sigma\\ shifts the disk (\\\sigma = 0\\: balanced, \\\sigma =
1\\: gain increases only up to the classical multiplicative form). It
equals \\\alpha\_{\max} = 1 / \\S + \frac{\sigma - 1}{2}\\\_\infty\\ and
guarantees *simultaneous* gain and phase variations inside the disk, a
stronger statement than the classical margins. The equivalent gain
margin interval and phase margin are reported. The nominal closed loop
must be stable.

## Usage

``` r
ctrl_disk_margin(L, skew = 0)
```

## Arguments

- L:

  SISO loop transfer function.

- skew:

  Skew \\\sigma\\.

## Value

A list of class `ctrl_disk_margin` with `alpha`, `gain_margin` (lower
and upper factors), `gain_margin_db`, `phase_margin` (degrees) and
`frequency` (the critical frequency).

## References

Seiler, P., Packard, A. and Gahinet, P. (2020). An introduction to disk
margins. *IEEE Control Systems Magazine*, 40(5), 78-95.

## Examples

``` r
L <- ctrl_tf(25, c(1, 10, 10, 10))
ctrl_disk_margin(L)
#> Disk margin (skew 0): alpha = 0.4581 at w = 1.955
#>   gain margin: [0.6273, 1.594] (-4.051 to 4.051 dB)
#>   phase margin: +/- 25.8 deg
```
