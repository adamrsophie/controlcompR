# Invariant zeros

For a square state-space model the invariant zeros are the finite
generalized eigenvalues of the Rosenbrock system pencil
\$\$\begin{pmatrix} A & B \\ C & D \end{pmatrix} - \lambda
\begin{pmatrix} I & 0 \\ 0 & 0 \end{pmatrix},\$\$ computed with the QZ
algorithm. A SISO transfer function returns its numerator roots. MIMO
transfer functions are converted to state space, so they must be square
as well; nonsquare models are not supported in this version.

## Usage

``` r
ctrl_zeros(sys)
```

## Arguments

- sys:

  A model.

## Value

A complex vector of zeros.

## Examples

``` r
ctrl_zeros(ctrl_tf(c(1, 2), c(1, 3, 2)))
#> [1] -2+0i
```
