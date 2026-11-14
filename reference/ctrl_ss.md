# State-space models

Construct a linear time-invariant state-space model \$\$\dot x(t) = A
x(t) + B u(t), \quad y(t) = C x(t) + D u(t)\$\$ in continuous time, or
\$\$x\_{k+1} = A x_k + B u_k, \quad y_k = C x_k + D u_k\$\$ in discrete
time with sample interval `dt`.

## Usage

``` r
ctrl_ss(
  A,
  B = NULL,
  C = NULL,
  D = NULL,
  dt = NULL,
  state_names = NULL,
  input_names = NULL,
  output_names = NULL
)
```

## Arguments

- A:

  State matrix (`n x n`).

- B:

  Input matrix (`n x m`).

- C:

  Output matrix (`p x n`).

- D:

  Feedthrough matrix (`p x m`); defaults to zero.

- dt:

  `NULL` or `0` for continuous time, or a positive sample interval.

- state_names:

  Optional unique state names.

- input_names:

  Optional unique input names.

- output_names:

  Optional unique output names. Defaults are `x1, x2, ...`, `u1, ...`
  and `y1, ...`.

## Value

An object of class `c("ctrl_ss", "ctrl_system")`.

## Details

The time domain is explicit: `dt = NULL` or `dt = 0` gives a continuous
time model, and a positive `dt` gives a discrete time model sampled
every `dt` time units. There is no "unknown sample time"; discrete
models always carry their sample interval.

Matrices are validated for finiteness and for consistent dimensions: `A`
is `n x n`, `B` is `n x m`, `C` is `p x n` and `D` is `p x m`. `D`
defaults to a zero matrix. A model with `n = 0` is a static gain
`y = D u`; pass `A = NULL` together with `D`, or use
[`ctrl_gain()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_gain.md).
Scalars are accepted for `1 x 1` matrices and vectors for single-row or
single-column matrices only when the shape is unambiguous.

## See also

[`ctrl_tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf.md),
[`ctrl_ssdata()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ssdata.md),
[`ctrl_c2d()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_c2d.md),
[`ctrl_poles()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_poles.md)

## Examples

``` r
A <- matrix(c(0, 1, -2, -3), 2, byrow = TRUE)
plant <- ctrl_ss(A, B = c(0, 1), C = c(1, 0))
plant
#> State-space model (continuous time)
#> 2 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1 x2
#> x1  0  1
#> x2 -2 -3
#> 
#> B =
#>    u1
#> x1  0
#> x2  1
#> 
#> C =
#>    x1 x2
#> y1  1  0
#> 
#> D =
#>    u1
#> y1  0
ctrl_poles(plant)
#> [1] -1+0i -2+0i

# A discrete time model sampled every 0.1 time units
ctrl_ss(0.9, 0.1, 1, dt = 0.1)
#> State-space model (discrete time, dt = 0.1)
#> 1 states, 1 inputs, 1 outputs
#> 
#> A =
#>     x1
#> x1 0.9
#> 
#> B =
#>     u1
#> x1 0.1
#> 
#> C =
#>    x1
#> y1  1
#> 
#> D =
#>    u1
#> y1  0
```
