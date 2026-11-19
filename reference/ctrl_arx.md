# ARX model estimation

Estimate the single-output ARX model \$\$y_k + a_1 y\_{k-1} + \cdots +
a\_{n_a} y\_{k-n_a} = \sum\_{i=1}^{m} \sum\_{j=0}^{n_b - 1} b\_{ij}
u\_{i, k - n_k - j} + e_k\$\$ by linear least squares (QR). The
convention for `nb` and `nk` matches the `sysid` package: `nb`
coefficients starting at delay `nk`. Standard errors assume white
equation errors; they are not valid when the true noise is colored (use
an output-error or Box-Jenkins model from `sysid` through
[`ctrl_from_sysid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_from_sysid.md)
in that case).

## Usage

``` r
ctrl_arx(data, na, nb, nk = 1, start = NULL)
```

## Arguments

- data:

  A `ctrl_iddata` object with one output.

- na:

  Number of output (denominator) coefficients.

- nb:

  Number of input coefficients per input.

- nk:

  Input delay in samples.

- start:

  First output sample used as a regression response. The default is the
  earliest sample for which all regressors exist. A later common `start`
  lets models of different orders be fitted and compared on the same
  responses (as
  [`ctrl_arx_order()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_arx_order.md)
  does).

## Value

An object of class `c("ctrl_arx", "ctrl_idmodel")` with `sys` (a
discrete `ctrl_tf` from inputs to the output), `a`, `b`, `vcov`,
`noise_var`, `residuals`, information criteria, `nobs` and `rows` (the
response samples used).

## Examples

``` r
set.seed(1)
u <- rnorm(300)
y <- as.numeric(stats::filter(0.5 * c(0, u[-300]), 0.8, "recursive")) +
  rnorm(300, sd = 0.05)
m <- ctrl_arx(ctrl_iddata(y, u), na = 1, nb = 1, nk = 1)
m$a
#> [1] -0.7985921
m$b
#>           [,1]
#> [1,] 0.4962837
```
