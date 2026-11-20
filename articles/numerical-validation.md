# Numerical validation, conventions and scope

``` r

library(controlcompR)
```

A control package is only useful if its numbers can be trusted. This
article states the conventions, explains how results are verified, lists
what is implemented, experimental, planned or unsupported, and describes
the limits of every guarantee.

## Conventions

| quantity | convention |
|----|----|
| continuous Lyapunov | A X + X A^\top + Q = 0 (same as python-control and MATLAB `lyap`) |
| discrete Lyapunov | A X A^\top - X + Q = 0 |
| CARE | A^\top P + P A - (PB + N) R^{-1} (B^\top P + N^\top) + Q = 0, stabilizing solution |
| DARE | A^\top P A - P - (A^\top P B + N)(R + B^\top P B)^{-1}(B^\top P A + N^\top) + Q = 0 |
| LQR law | u = -Kx |
| Kalman `P` | prior covariance P\_{k\|k-1}; `Z` is the posterior P\_{k\|k} |
| Kalman `M`, `L` | measurement update gain and one-step predictor gain |
| feedback | `ctrl_feedback(G, K)` is (I + GK)^{-1}G, negative by default |
| series | `ctrl_series(G1, G2)` is G_2 G_1 (signal flows from 1 to 2) |
| margins | loop gain L under unity negative feedback |
| frequency units | rad/s (rad/sample divided by `dt` for discrete models) |
| discrete time | `dt > 0` is the sample interval; `dt = 0` is continuous |
| MPC decision vector | (x_0, \ldots, x_N, u_0, \ldots, u\_{N-1}), u_k acts between x_k and x\_{k+1} |
| H-infinity norm | peak of the largest singular value over frequency; `Inf` for unstable models |
| H2 norm | `Inf` for unstable models or continuous models with direct feedthrough |

## Algorithms

Numerically sensitive computations use orthogonal decompositions and
avoid explicitly forming squared or ill-conditioned quantities. Many
results are checked after the fact, but the checks detect failures; they
do not prove backward stability for every input. The main methods:

- **Riccati equations**: ordered real Schur decomposition of the
  Hamiltonian (continuous) or symplectic pencil (discrete, generalized
  Schur), with the stabilizing subspace selected by eigenvalue location.
  The relative residual (accepted up to 1e-6) and the closed-loop
  eigenvalues are verified; a solution that fails either check raises an
  error instead of being returned.
- **Lyapunov and Stein equations**: real Schur form, then
  Bartels-Stewart (continuous) or block back substitution (discrete).
- **Controllability and observability**: orthogonal staircase forms
  (never the rank of the controllability matrix, which is numerically
  meaningless beyond small systems).
- **Gramians and Hankel singular values**: Cholesky factors of the
  Gramians computed directly by Hammarling’s method on the complex Schur
  form, so the Gramians are never formed and their condition number is
  never squared; the H2 norm is \\C L_c\\\_F with the same factor.
- **H-infinity norm**: the Bruinsma-Steinbuch level-set iteration on the
  Hamiltonian, with the imaginary-axis eigenvalue test.
- **Discretization**: matrix exponential of an augmented matrix
  (zero-order and first-order hold), bilinear (Tustin) and matched
  pole-zero transformations.
- **Frequency response**: a complex linear solve per frequency on a
  diagonally balanced state-space model. Transfer function objects are
  evaluated from their polynomial coefficients, which is accurate for
  low orders; convert high-order transfer functions to state space first
  ([`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md))
  for better conditioning.
- **Pole placement**: robust eigenstructure assignment (Kautsky, Nichols
  and Van Dooren), with the eigenvalues of the result verified.
- **QPs**: OSQP with tight tolerances and solution polishing; every MPC
  step reports the solver status, bound violation and dynamics residual.

``` r

A <- matrix(c(0, 1, -2, -3), 2, byrow = TRUE)
B <- matrix(c(0, 1))
sol <- ctrl_care(A, B, diag(2), 1)
sol$residual
#> [1] 1.720524e-15
sol$poles
#> [1] -1.000000+0i -2.236068+0i
```

## Tolerances

Rank decisions and structural tests use relative tolerances based on the
size and the norm of the problem, of the form \mathrm{tol} = c\\ n^2\\
\varepsilon\\ \\M\\ for the staircase forms. Stability is decided on a
state realization balanced together with its input and output matrices,
with per-pole bands from the first-order eigenvalue error bound, so the
decision does not depend on state units. Covariance matrices are
inverted on their support after equilibration (scaling to unit
diagonal), so a small but positive variance is never treated as zero.
Boundary cases are reported (`"marginally stable"` for linear models,
`"inconclusive"` for equilibria of nonlinear models) instead of being
forced into stable or unstable. Every function that takes a tolerance
documents it.

## Testing strategy

The test suite checks results against properties that must hold, rather
than against stored numbers:

- residuals of matrix equations and the eigenvalues of closed loops;
- invariance under similarity transformations (poles, zeros, Hankel
  singular values, norms);
- identities: series and feedback algebra, the separation principle for
  LQG, `ctrl_ss2tf(ctrl_tf2ss(G))` equals `G`;
- exact cross-checks between independent algorithms: MPC without active
  constraints equals LQR with the Riccati terminal weight, MHE without
  active constraints equals the Kalman filter, the H2 norm from the
  Gramian factor equals the integral of the frequency response, the
  H-infinity norm equals the peak of a dense frequency sweep;
- analytic results: first and second order step responses,
  Ziegler-Nichols tables, balanced truncation error bounds;
- invariance under changes of units: diagonal and orthogonal state
  scalings, input and cost scalings and time scalings leave Riccati
  solutions, Hankel singular values, norms, stability decisions, Kalman,
  UKF, smoother and MHE posteriors unchanged (including covariances
  whose entries span many orders of magnitude);
- rejection of invalid problems: indefinite weights, unstabilizable
  pairs, singular measurement noise, ill-posed interconnections,
  infeasible QPs, nonfinite inputs and arithmetic that overflows.

``` r

G <- ctrl_ss(matrix(c(-1, 2, 0, -3), 2), matrix(c(1, 0)), matrix(c(0, 1), 1))
T <- matrix(c(2, 1, 1, 1), 2)
G2 <- ctrl_ss(T %*% G$A %*% solve(T), T %*% G$B, G$C %*% solve(T))
c(ctrl_norm(G, "inf"), ctrl_norm(G2, "inf"))
#> [1] 0.6666667 0.6666667
rbind(ctrl_hsvd(G), ctrl_hsvd(G2))
#>           [,1]       [,2]
#> [1,] 0.3871459 0.05381261
#> [2,] 0.3871459 0.05381261
```

## Reference comparison

The repository contains `benchmarks/reference/`: a set of problems
(`problems.json`), a Python script that solves them with python-control
and slycot (the SLICOT Fortran library), and an R script that compares
controlcompR’s results. The comparison runs on GitHub Actions; it covers
continuous and discrete Riccati and Lyapunov equations, H-infinity and
H2 norms, Hankel singular values, zero-order hold discretization, pole
placement and H-infinity synthesis. Relative errors are required to be
below 10^{-7} (and 10^{-3} for the bisection level of H-infinity
synthesis, where both packages stop at a tolerance).

## Capability matrix

| area | status | notes |
|----|----|----|
| state space, transfer function, zero-pole-gain models | implemented | MIMO; conversions verified by round trips |
| discretization (zoh, foh, tustin, matched) | implemented |  |
| interconnections and block diagrams | implemented | algebraic loops checked for well-posedness |
| poles, zeros, controllability, observability, minimal realization | implemented | staircase forms |
| Lyapunov and Riccati equations | implemented | residual and stability checks |
| time responses (step, impulse, initial, arbitrary input) | implemented | exact zero-order hold for sampled inputs |
| frequency responses, margins, disk margins, norms | implemented |  |
| LQR (continuous, discrete, finite horizon) | implemented |  |
| pole placement and observers | implemented |  |
| Kalman filter, smoother, LQG | implemented |  |
| PID with filter, anti-windup, bumpless transfer, tuning rules | implemented |  |
| linear MPC (bounds, rates, soft constraints, terminal sets) | implemented | OSQP |
| offset-free MPC and MHE | implemented |  |
| balanced truncation and residualization | implemented | stable part; unstable part kept |
| ARX, subspace identification, validation tools | implemented | other model structures through sysid |
| H2 synthesis | implemented | continuous time, D\_{11} = 0 |
| H-infinity synthesis | experimental | continuous time, D\_{11} = 0, every result checked against the upper bound of an independent norm computation |
| nonlinear simulation, linearization, equilibria | implemented |  |
| EKF and UKF | implemented |  |
| nonlinear MPC | experimental | single shooting, local optimum |
| discrete time H-infinity synthesis | planned | refused with an error |
| general D\_{11} in H-infinity synthesis | planned | refused with an error |
| structured singular value, mu synthesis | unsupported |  |
| explicit MPC, robust tube MPC | unsupported |  |
| descriptor (E \ne I) systems | unsupported |  |
| hybrid and switched systems | unsupported |  |

## What the guarantees cover

A guarantee in this package is a mathematical statement about a model.
LQR optimality, MPC constraint satisfaction, H-infinity bounds and
robust stability tests all hold for the model and the uncertainty
description that were given. None of them covers unmodeled dynamics,
sensor faults, timing jitter, quantization or a wrong model. Verify
every design in simulation with realistic disturbances and noise, and
treat deployment to a physical system as an engineering task with its
own safety analysis. This package does not generate certified code for
embedded targets.
