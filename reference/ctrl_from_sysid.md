# Convert a `sysid` model

Convert a model estimated by the `sysid` package (an `estpoly` fit or an
`idpoly` object) to `controlcompR` transfer functions. The `sysid`
polynomial model \$\$A(q) y = \frac{B(q)}{F(q)} q^{-n_k} u +
\frac{C(q)}{D(q)} e\$\$ (polynomials in \\q^{-1}\\) gives the plant
\\G(z) = z^{-n_k} B(z^{-1}) / (A(z^{-1}) F(z^{-1}))\\ and the noise
model \\H(z) = C(z^{-1}) / (A(z^{-1}) D(z^{-1}))\\.

## Usage

``` r
ctrl_from_sysid(model)
```

## Arguments

- model:

  An `estpoly` or `idpoly` object from `sysid`.

## Value

A list with `plant` and `noise` (`ctrl_tf`), `noise_var` and `type`.

## Examples

``` r
if (requireNamespace("sysid", quietly = TRUE)) {
  set.seed(1)
  u <- rnorm(200)
  y <- as.numeric(stats::filter(0.5 * c(0, u[-200]), 0.8, "recursive")) +
    rnorm(200, sd = 0.05)
  fit <- sysid::arx(sysid::idframe(output = y, input = u, Ts = 1), c(1, 1, 1))
  ctrl_from_sysid(fit)$plant
}
#> Transfer function (discrete time, dt = 1)
#> 
#>     0.4977
#>   ----------
#>   z - 0.7979
```
