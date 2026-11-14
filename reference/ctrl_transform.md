# Similarity transformation of a state-space model

Apply the change of state coordinates \\z = T x\\, giving \\(T A T^{-1},
T B, C T^{-1}, D)\\.

## Usage

``` r
ctrl_transform(sys, T)
```

## Arguments

- sys:

  A `ctrl_ss` model.

- T:

  Nonsingular `n x n` transformation matrix.

## Value

The transformed `ctrl_ss` model (input-output behavior unchanged).

## Examples

``` r
sys <- ctrl_ss(matrix(c(0, 1, -2, -3), 2, byrow = TRUE), c(0, 1), c(1, 0))
ctrl_transform(sys, matrix(c(1, 1, 0, 1), 2))
#> State-space model (continuous time)
#> 2 states, 1 inputs, 1 outputs
#> 
#> A =
#>    x1 x2
#> x1 -1  1
#> x2  0 -2
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
```
