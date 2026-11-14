# Connect models by signal names

Build an interconnection of several models whose input and output names
define the wiring, as in a block diagram. Every block input whose name
equals a block output name is driven by that output; the remaining block
inputs must be listed in `inputs` and become external inputs. Use
[`ctrl_sumblock()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_sumblock.md)
for sums of signals. Output names must be unique across blocks.

## Usage

``` r
ctrl_connect(..., inputs, outputs)
```

## Arguments

- ...:

  Models with named signals.

- inputs:

  Names of the external inputs.

- outputs:

  Names of the block outputs to expose.

## Value

A `ctrl_ss` model from `inputs` to `outputs`.

## Details

Algebraic loops through feedthrough terms are resolved exactly; an error
is raised when the loop matrix is singular (ill-posed diagram).

## Examples

``` r
G <- ctrl_tf(1, c(1, 1), input_names = "u", output_names = "y")
K <- ctrl_tf(5, 1, input_names = "e", output_names = "u")
S <- ctrl_sumblock("e = r - y")
T <- ctrl_connect(G, K, S, inputs = "r", outputs = "y")
ctrl_dcgain(T)  # 5 / 6
#>           r
#> y 0.8333333
```
