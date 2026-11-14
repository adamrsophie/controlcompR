# Models, conversions and interconnections

``` r

library(controlcompR)
```

## State-space models

A state-space model is \dot x = A x + B u, \qquad y = C x + D u in
continuous time or x\_{k+1} = A x_k + B u_k, \qquad y_k = C x_k + D u_k
in discrete time with sample interval `dt`.
[`ctrl_ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ss.md)
validates the dimensions, requires finite real matrices and accepts
vectors for single inputs and outputs:

``` r

sys <- ctrl_ss(matrix(c(0, 1, -4, -0.8), 2, byrow = TRUE), B = c(0, 1),
               C = c(1, 0))
ctrl_dims(sys)
#>  states  inputs outputs 
#>       2       1       1
ctrl_is_continuous(sys)
#> [1] TRUE
sysd <- ctrl_ss(0.9, 0.1, 1, dt = 0.5)
ctrl_sample_time(sysd)
#> [1] 0.5
```

`D` defaults to zero. Signal names are carried through conversions and
interconnections:

``` r

sys2 <- ctrl_ss(diag(c(-1, -2)), diag(2), diag(2),
                input_names = c("heater", "fan"),
                output_names = c("temperature", "humidity"))
ctrl_signal_names(sys2)
#> $states
#> [1] "x1" "x2"
#> 
#> $inputs
#> [1] "heater" "fan"   
#> 
#> $outputs
#> [1] "temperature" "humidity"
sys2[1, 2]  # output 1, input 2
#> State-space model (continuous time)
#> 2 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1 x2
#> x1 -1  0
#> x2  0 -2
#> 
#> B =
#>    fan
#> x1   0
#> x2   1
#> 
#> C =
#>             x1 x2
#> temperature  1  0
#> 
#> D =
#>             fan
#> temperature   0
```

A static gain is a model without states:

``` r

ctrl_gain(matrix(c(1, 2, 3, 4), 2))
#> State-space model (continuous time)
#> 0 states, 2 inputs, 2 outputs
#> 
#> A =
#>   (empty 0 x 0)
#> 
#> B =
#>   (empty 0 x 2)
#> 
#> C =
#>   (empty 2 x 0)
#> 
#> D =
#>    u1 u2
#> y1  1  3
#> y2  2  4
```

## Transfer functions and zero-pole-gain models

Polynomial coefficients are in descending powers of `s` (or `z`):

``` r

G <- ctrl_tf(c(1, 3), c(1, 2, 5))
G
#> Transfer function (continuous time)
#> 
#>         s +     3
#>   ---------------------
#>   s^2 +     2 s +     5
ctrl_zeros(G)
#> [1] -3+0i
ctrl_poles(G)
#> [1] -1+2i -1-2i
```

Multivariable transfer functions are list matrices of numerators and
denominators. The zero-pole-gain form is often better conditioned for
high order models:

``` r

H <- ctrl_zpk(zeros = -3, poles = c(-1 + 2i, -1 - 2i), gain = 1)
H
#> Zero-pole-gain model (continuous time)
#>   gain: 1
#>   zeros: -3+0i
#>   poles: -1+2i, -1-2i
ctrl_zpk2tf(H)
#> Transfer function (continuous time)
#> 
#>         s +     3
#>   ---------------------
#>   s^2 +     2 s +     5
```

## Conversions

[`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
builds a realization (minimal by default) and
[`ctrl_ss2tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
computes the transfer function of each channel:

``` r

ss <- ctrl_tf2ss(G)
ss
#> State-space model (continuous time)
#> 2 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1 x2
#> x1 -2 -5
#> x2  1  0
#> 
#> B =
#>    u1
#> x1  1
#> x2  0
#> 
#> C =
#>    x1 x2
#> y1  1  3
#> 
#> D =
#>    u1
#> y1  0
ctrl_ss2tf(ss)
#> Transfer function (continuous time)
#> 
#>          1 s +     3
#>   ---------------------
#>   s^2 +     2 s +     5
```

Realizations are not unique: any similarity transformation T gives the
same input-output behavior, and
[`ctrl_transform()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_transform.md)
applies one explicitly. Pole-zero cancellations are removed by
[`ctrl_minreal()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_minreal.md),
whose default tolerance is scaled by the norms of the model:

``` r

Gc <- ctrl_tf(c(1, 1), c(1, 3, 2))  # (s + 1) / ((s + 1)(s + 2))
ctrl_minreal(ctrl_tf2ss(Gc, minimal = FALSE))
#> State-space model (continuous time)
#> 1 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1
#> x1 -2
#> 
#> B =
#>        u1
#> x1 0.7071
#> 
#> C =
#>       x1
#> y1 1.414
#> 
#> D =
#>    u1
#> y1  0
```

## Discretization

[`ctrl_c2d()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_c2d.md)
supports exact zero-order hold (`"zoh"`, the default), first-order hold
(`"foh"`), the bilinear transform (`"tustin"`, with optional frequency
prewarping) and the matched pole-zero method. The zero-order hold is
computed from the exponential of an augmented matrix \exp\begin{pmatrix}
A & B \\ 0 & 0 \end{pmatrix} T, which avoids inverting `A`:

``` r

Gd <- ctrl_c2d(sys, 0.1)
Gd
#> State-space model (discrete time, dt = 0.1)
#> 2 states, 1 inputs, 1 outputs
#> 
#> A =
#>         x1      x2
#> x1  0.9806 0.09547
#> x2 -0.3819 0.90420
#> 
#> B =
#>          u1
#> x1 0.004853
#> x2 0.095470
#> 
#> C =
#>    x1 x2
#> y1  1  0
#> 
#> D =
#>    u1
#> y1  0
ctrl_poles(Gd)
#> [1] 0.9424012+0.1870729i 0.9424012-0.1870729i
exp(ctrl_poles(sys) * 0.1)  # zoh maps poles by exp(sT)
#> [1] 0.9424012+0.1870729i 0.9424012-0.1870729i
```

The discretizations differ in what they preserve. The zero-order hold is
exact for piecewise constant inputs; Tustin preserves stability and maps
the whole imaginary axis onto the unit circle, with frequency warping
corrected at one frequency by `prewarp`:

``` r

w0 <- 2
tus <- ctrl_c2d(sys, 0.1, method = "tustin", prewarp = w0)
c(continuous = Mod(ctrl_freqresp(sys, w0)$response[1, 1, 1]),
  tustin = Mod(ctrl_freqresp(tus, w0)$response[1, 1, 1]))
#> continuous     tustin 
#>      0.625      0.625
```

[`ctrl_d2c()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_d2c.md)
inverts the zero-order hold and Tustin maps (the zero-order hold inverse
needs a matrix logarithm and fails for poles on the negative real axis).

## Interconnections

Series, parallel and feedback connections work on any mix of models in
the same time domain. Arithmetic operators are shorthand:

``` r

C <- ctrl_tf(c(2, 1), c(1, 0))   # PI controller
L <- C * G                        # series: G after C
T <- ctrl_feedback(L)             # unity negative feedback
ctrl_poles(T)
#> [1] -0.2731763+0.000000i -1.8634119+2.740368i -1.8634119-2.740368i
S <- ctrl_feedback(1, L)          # sensitivity 1 / (1 + L)
```

`ctrl_series(sys1, sys2)` applies `sys1` first, so `C * G` equals
`ctrl_series(G, C)` for single-input single-output models. For block
diagrams with named signals,
[`ctrl_sumblock()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_sumblock.md)
and
[`ctrl_connect()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_connect.md)
wire blocks by signal name:

``` r

Gs <- ctrl_tf2ss(G)
Gs$input_names <- "u"
Gs$output_names <- "y"
Cs <- ctrl_tf2ss(C)
Cs$input_names <- "e"
Cs$output_names <- "u"
sum_e <- ctrl_sumblock("e = r - y")
cl <- ctrl_connect(Gs, Cs, sum_e, inputs = "r", outputs = c("y", "u"))
cl
#> State-space model (continuous time)
#> 3 states, 1 inputs, 2 outputs
#> 
#> A =
#>    x1  x2 x3
#> x1 -4 -11  1
#> x2  1   0  0
#> x3 -1  -3  0
#> 
#> B =
#>    r
#> x1 2
#> x2 0
#> x3 1
#> 
#> C =
#>   x1 x2 x3
#> y  1  3  0
#> u -2 -6  1
#> 
#> D =
#>   r
#> y 0
#> u 2
```

The result is the closed loop from `r` to `y` and `u`, equal to the
feedback connection above:

``` r

ctrl_dcgain(cl)
#>          r
#> y 1.000000
#> u 1.666667
```
