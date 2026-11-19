# Residual analysis for identified models

Compute the residuals of a model on data (one-step prediction errors for
ARX models, output errors otherwise) and test them: the autocorrelation
should be within the 99 percent band \\\pm 2.58/\sqrt{N}\\ for white
residuals (Ljung-Box test reported), and the cross-correlation with the
input should vanish at nonnegative lags (residuals independent of
current and past inputs; `ccf` reports \\\mathrm{cor}(e\_{t+k}, u_t)\\).
All outputs and input-output pairs are tested. Output-error models are
not expected to have white residuals; the cross-correlation test is the
relevant one for them.

## Usage

``` r
ctrl_resid(model, data, lags = 20L)
```

## Arguments

- model:

  An identified model.

- data:

  Data (preferably validation data).

- lags:

  Maximum lag; capped at the number of residuals minus one.

## Value

A list of class `ctrl_resid` with `residuals`, `acf`, `ccf` (data frames
with output, input, lag, value and band), `ljung_box` (p-value per
output; `NA` when the lags do not exceed the fitted parameters), `lags`
(the lags used) and flags `white` and `independent`.

## Examples

``` r
set.seed(1)
u <- rnorm(300)
y <- as.numeric(stats::filter(0.5 * c(0, u[-300]), 0.8, "recursive")) +
  rnorm(300, sd = 0.05)
d <- ctrl_iddata(y, u)
r <- ctrl_resid(ctrl_arx(d, 1, 1, 1), d)
r$white
#> [1] FALSE
```
