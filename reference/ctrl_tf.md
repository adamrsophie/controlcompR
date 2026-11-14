# Transfer-function models

Construct a rational transfer function \$\$G(s) = \frac{b_m s^m +
\cdots + b_1 s + b_0}{a_n s^n + \cdots + a_1 s + a_0}\$\$ in continuous
time (`dt = NULL` or `0`) or in the variable `z` for a discrete time
model with sample interval `dt`.

## Usage

``` r
ctrl_tf(num, den = 1, dt = NULL, input_names = NULL, output_names = NULL)
```

## Arguments

- num:

  Numerator coefficients (descending powers), a list matrix for MIMO
  models, or the string `"s"` or `"z"` for the Laplace or shift
  variable.

- den:

  Denominator coefficients (descending powers) or a list matrix.

- dt:

  `NULL` or `0` for continuous time, or a positive sample interval.

- input_names:

  Optional unique input names.

- output_names:

  Optional unique output names. Defaults are `x1, x2, ...`, `u1, ...`
  and `y1, ...`.

## Value

An object of class `c("ctrl_tf", "ctrl_system")` with fields `num` and
`den` (list matrices), `dt` and signal names.

## Details

**Coefficient ordering.** Polynomials are given in *descending* powers,
as in MATLAB and python-control: `c(1, 3, 2)` is \\s^2 + 3s + 2\\.
Leading zeros are removed and the denominator is normalized to be monic;
the numerator is scaled by the same factor.

**MIMO models.** For `p` outputs and `m` inputs, pass `num` as a `p x m`
list matrix of numerator vectors and `den` either as a single vector
(common denominator) or as a `p x m` list matrix.

**Improper models.** A transfer function with a numerator degree above
the denominator degree can be created (for example `ctrl_tf("s")`, which
is handy for building models with arithmetic) and evaluated in the
frequency domain, but it has no finite-dimensional state-space
realization and
[`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
refuses to convert it.

## See also

[`ctrl_zpk()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_zpk.md),
[`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md),
[`ctrl_ss2tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)

## Examples

``` r
G <- ctrl_tf(1, c(1, 3, 2))
G
#> Transfer function (continuous time)
#> 
#>               1
#>   ---------------------
#>   s^2 +     3 s +     2
s <- ctrl_tf("s")
H <- (s + 1) / (s^2 + 2 * s + 5)
H
#> Transfer function (continuous time)
#> 
#>         s +     1
#>   ---------------------
#>   s^2 +     2 s +     5
# A discrete time integrator with sample interval 0.1
ctrl_tf(c(0.1, 0), c(1, -1), dt = 0.1)
#> Transfer function (discrete time, dt = 0.1)
#> 
#>      0.1 z
#>   ---------
#>   z -     1
```
