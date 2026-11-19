# ARX model order selection

Fit ARX models over a grid of orders and report information criteria on
the estimation data and, when given, simulation fit on validation data.
Prefer validation fit and residual tests to a single criterion; AIC
tends to overfit, BIC is more conservative.

## Usage

``` r
ctrl_arx_order(data, na = 1:4, nb = 1:4, nk = 1, validation = NULL)
```

## Arguments

- data:

  Estimation data.

- na, nb:

  Candidate orders.

- nk:

  Delay.

- validation:

  Optional validation data.

## Value

A data frame sorted by BIC, with `nobs`. All candidates are fitted on
the same response samples, so their criteria are comparable.

## Examples

``` r
set.seed(1)
u <- rnorm(400)
y <- as.numeric(stats::filter(0.5 * c(0, u[-400]), 0.8, "recursive")) +
  rnorm(400, sd = 0.05)
parts <- ctrl_idsplit(ctrl_iddata(y, u))
head(ctrl_arx_order(parts$estimation, 1:3, 1:3, validation = parts$validation))
#>   na nb nk npar nobs       aic       bic validation_fit
#> 9  3  3  1    6  277 -778.8647 -757.1206       92.71684
#> 5  2  2  1    4  277 -755.7826 -741.2866       92.72494
#> 6  3  2  1    5  277 -754.2925 -736.1724       92.71098
#> 8  2  3  1    5  277 -753.7902 -735.6701       92.72485
#> 1  1  1  1    2  277 -681.9408 -674.6927       92.54913
#> 2  2  1  1    3  277 -684.9317 -674.0596       92.64053
```
