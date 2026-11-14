# Summing junction

Create a static gain block from an equation such as `"e = r - y"`. Use
it with
[`ctrl_connect()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_connect.md)
to build block diagrams by signal name.

## Usage

``` r
ctrl_sumblock(expr, dt = NULL)
```

## Arguments

- expr:

  A string `"out = +a - b + c"`.

- dt:

  Time domain (a static block adopts the domain of the others).

## Value

A `ctrl_ss` static gain with named inputs and output.

## Examples

``` r
ctrl_sumblock("e = r - y")
#> State-space model (continuous time)
#> 0 states, 2 inputs, 1 outputs
#> 
#> A =
#>   (empty 0 x 0)
#> 
#> B =
#>   (empty 0 x 2)
#> 
#> C =
#>   (empty 1 x 0)
#> 
#> D =
#>   r  y
#> e 1 -1
```
