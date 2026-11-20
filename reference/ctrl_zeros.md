# Invariant zeros

For a square state-space model the invariant zeros are the finite
generalized eigenvalues of the Rosenbrock system pencil
\$\$\begin{pmatrix} A & B \\ C & D \end{pmatrix} - \lambda
\begin{pmatrix} I & 0 \\ 0 & 0 \end{pmatrix},\$\$ computed with the QZ
algorithm. When `D` is well-conditioned relative to its own size
(reciprocal condition number above \\\sqrt\epsilon\\), there are exactly
`n` finite zeros; zeros too large for the pencil to separate from
infinity (from a tiny feedthrough) are then taken as the largest
eigenvalues of \\A - B D^{-1} C\\, so they are kept. Otherwise a
generalized eigenvalue that cannot be told apart from infinity at
working precision is treated as infinite with a warning. Zeros need a
regular pencil (full normal rank); this is checked by evaluating the
transfer matrix at points away from every computed eigenvalue and pole,
and a singular pencil raises an error. A SISO transfer function returns
its numerator roots. MIMO transfer functions are converted to state
space, so they must be square as well; nonsquare models are not
supported in this version.

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
