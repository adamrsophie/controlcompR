# Split identification data into estimation and validation parts

Split identification data into estimation and validation parts

## Usage

``` r
ctrl_idsplit(data, fraction = 0.7)
```

## Arguments

- data:

  A `ctrl_iddata` object.

- fraction:

  Fraction of samples used for estimation (the first part; the remainder
  is kept for validation, preserving time order).

## Value

A list with `estimation` and `validation` data sets.

## Examples

``` r
d <- ctrl_iddata(rnorm(100), rnorm(100))
parts <- ctrl_idsplit(d, 0.7)
```
