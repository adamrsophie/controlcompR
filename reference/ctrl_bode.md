# Bode, Nyquist and Nichols data

`ctrl_bode()` returns magnitude (absolute and dB) and phase (degrees)
over a frequency grid. Phase is unwrapped along the grid: the result is
continuous only when the grid resolves every rapid phase change, so
check lightly damped modes with a denser `w`. `ctrl_nyquist()` returns
the complex response for positive and negative frequencies and, for SISO
models without poles on the stability boundary, the winding number of
\\1 + L\\ around the origin (the number of clockwise encirclements of
\\-1\\); with `P` open-loop unstable poles the closed loop has \\Z = N +
P\\ unstable poles. The winding number is computed on its own contour:
for continuous models the imaginary axis is extended until the remaining
tail of \\L - D\\ cannot reach \\-1\\, and the contour is refined until
the phase of \\1 + L\\ changes by less than \\\pi/8\\ between samples.
When this cannot be established (for example \\1 + D = 0\\),
`encirclements` is `NA`. `ctrl_nichols()` returns phase against
magnitude in dB.

## Usage

``` r
ctrl_bode(sys, w = NULL)

ctrl_nyquist(sys, w = NULL)

ctrl_nichols(sys, w = NULL)
```

## Arguments

- sys:

  A model.

- w:

  Frequencies (rad per time unit).

## Value

Objects of class `ctrl_bode`, `ctrl_nyquist` or `ctrl_nichols` with a
`data` data frame and a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.

## Examples

``` r
L <- ctrl_tf(10, c(1, 3, 3, 1))
b <- ctrl_bode(L)
head(b$data)
#>            w magnitude magnitude_db     phase input output
#> 1 0.01000000  9.998500     19.99870 -1.718816    u1     y1
#> 2 0.01023352  9.998429     19.99864 -1.758951    u1     y1
#> 3 0.01047249  9.998355     19.99857 -1.800023    u1     y1
#> 4 0.01071705  9.998277     19.99850 -1.842054    u1     y1
#> 5 0.01096731  9.998196     19.99843 -1.885067    u1     y1
#> 6 0.01122342  9.998111     19.99836 -1.929083    u1     y1
ny <- ctrl_nyquist(L)
ny$encirclements
#> [1] 2
```
