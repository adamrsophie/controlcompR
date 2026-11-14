# Select inputs and outputs of a model

`sys[i, j]` keeps outputs `i` and inputs `j` (indices or names).

## Usage

``` r
# S3 method for class 'ctrl_ss'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A `ctrl_ss` model.

- i:

  Outputs to keep.

- j:

  Inputs to keep.

- ...:

  Unused.

- drop:

  Unused; dimensions are never dropped.

## Value

A `ctrl_ss` model.

## Examples

``` r
sys <- ctrl_ss(-diag(2), diag(2), diag(2))
sys[1, 2]
#> State-space model (continuous time)
#> 2 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1 x2
#> x1 -1  0
#> x2  0 -1
#> 
#> B =
#>    u2
#> x1  0
#> x2  1
#> 
#> C =
#>    x1 x2
#> y1  1  0
#> 
#> D =
#>    u2
#> y1  0
```
