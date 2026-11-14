# Compare identified models with data

Compute normalized fit percentages, \\100 (1 - \\y - \hat y\\ / \\y -
\bar y\\)\\, for the simulated output (input only, the demanding test
relevant for control design) and, for ARX models, the one-step-ahead
prediction. Use validation data that were not used for estimation. A
high prediction fit does not establish that a model is adequate for
feedback design.

## Usage

``` r
ctrl_compare(data, ...)
```

## Arguments

- data:

  A `ctrl_iddata` object (preferably validation data).

- ...:

  Identified models (`ctrl_idmodel` objects or plain models), optionally
  named.

## Value

A data frame with one row per model and output, with class
`ctrl_idcompare` (the simulated outputs are attached as attribute
`outputs` for plotting).

## Examples

``` r
set.seed(1)
u <- rnorm(400)
y <- as.numeric(stats::filter(0.5 * c(0, u[-400]), 0.8, "recursive")) +
  rnorm(400, sd = 0.05)
parts <- ctrl_idsplit(ctrl_iddata(y, u))
m1 <- ctrl_arx(parts$estimation, 1, 1, 1)
m2 <- ctrl_arx(parts$estimation, 2, 2, 1)
ctrl_compare(parts$validation, first = m1, second = m2)
#>    model output simulation_fit prediction_fit
#> 1  first     y1       92.55679       90.70529
#> 2 second     y1       92.72496       91.48786
```
