# controlcompR: Control Systems, Optimal Control, and Model Predictive Control in R

***controlcompR*** is an R framework for modeling, analyzing, designing
and simulating control systems. It follows the path

**model, analysis, design, estimation, constrained control,
validation**,

in the spirit of [python-control](https://python-control.org) and
[ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl),
but everything runs in R. Dense linear algebra uses LAPACK through a
small C layer (ordered real and generalized Schur forms), quadratic
programs are solved with [OSQP](https://osqp.org), and ODEs with
[deSolve](https://cran.r-project.org/package=deSolve).

Every numerically sensitive result is checked before it is returned:
Riccati solutions are verified by their residual and closed-loop
eigenvalues, H-infinity controllers are certified by an independent norm
computation, MPC steps report solver status, bound violations and
dynamics residuals, and problems the package cannot solve reliably are
refused with an error rather than approximated silently.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("adamrsophie/controlcompR")
```

## A first example

``` r

library(controlcompR)

# a mass with friction, force input, position output
G <- ctrl_ss(matrix(c(0, 1, 0, -0.5), 2, byrow = TRUE), c(0, 1), c(1, 0))

# LQR design, closed loop and step response
lqr <- ctrl_lqr(G, Q = diag(c(10, 1)), R = 0.1)
cl <- ctrl_closed_loop(G, lqr)
N <- ctrl_reference_gain(G, lqr)
plot(ctrl_step(ctrl_series(ctrl_gain(N), cl)))

# a PI loop on the transfer function: margins and peak sensitivity
P <- ctrl_tf(1, c(1, 0.5, 0))
C <- ctrl_tf(c(2, 0.2), c(1, 0))
ctrl_margin(C * P)
ctrl_norm(ctrl_feedback(1, C * P), "inf")
```

Constrained model predictive control of the discretized model:

``` r

Gd <- ctrl_c2d(G, 0.1)
mpc <- ctrl_mpc(Gd, horizon = 20, Q = diag(c(10, 1)), R = 0.1,
                u_min = -1, u_max = 1, x_min = c(-Inf, -0.5), x_max = c(Inf, 0.5))
sim <- ctrl_simulate_mpc(mpc, x0 = c(2, 0), steps = 80)
summary(sim)
plot(sim)
```

## What is implemented

| Area | Functionality |
|----|----|
| Models | state space, transfer functions, zero-pole-gain, static gains; MIMO; signal names; continuous and discrete time; conversions; zoh, foh, Tustin and matched discretization; `d2c` |
| Interconnections | series, parallel, feedback (algebraic loops checked), append, named-signal block diagrams, sum blocks, closed loops, sensitivity functions, linear fractional transformations |
| Analysis | poles, transmission zeros, damping, DC gain, stability with tolerances and marginal cases, staircase controllability and observability, Gramians, stabilizability, minimal realizations |
| Matrix equations | Lyapunov, Stein and Sylvester equations; continuous and discrete Riccati equations with cross terms and verified stabilizing solutions |
| Time responses | step, impulse, initial condition and arbitrary input responses with exact sampling; step metrics |
| Frequency responses | Bode, Nyquist, Nichols, singular values, gain, phase, delay and disk margins, bandwidth, H2 and H-infinity norms, root locus, pole-zero maps |
| State feedback | continuous, discrete and finite-horizon LQR, reference prefilters, robust pole placement |
| Estimation | Luenberger observers, steady-state and time-varying Kalman filters (Joseph form, missing data, NIS, likelihood), RTS smoother, LQG |
| PID | parallel and standard forms, derivative filter, setpoint weights, saturation, back-calculation and clamping anti-windup, bumpless transfer, Ziegler-Nichols, SIMC and AMIGO tuning, FOPDT fitting |
| MPC | sparse QP with input, rate, state, output and terminal constraints, soft constraints, control horizon, preview references, offset-free tracking, explicit failure handling |
| Estimation with constraints | moving horizon estimation with warm starts |
| Robust control | H2 synthesis, H-infinity synthesis (experimental, certified), mixed-sensitivity plants, small-gain robust stability and robust performance tests |
| Model reduction | Hankel singular values, balanced realization, balanced truncation and residualization with error bounds |
| Identification | ARX, order selection, subspace (MOESP) identification, simulation and prediction fit, residual tests, conversion of `sysid` models |
| Nonlinear systems | simulation, linearization, equilibria, discretization, EKF, UKF, nonlinear MPC (experimental), trajectory optimization |

Not available: discrete-time H-infinity synthesis, general `D11` in
H-infinity synthesis, structured singular values and mu synthesis,
explicit MPC, tube MPC, descriptor systems and hybrid systems. The
[numerical validation
article](https://adamrsophie.github.io/controlcompR/articles/numerical-validation.html)
lists the scope and conventions in detail.

## Validation

The test suite checks properties rather than stored numbers: residuals
of matrix equations, invariance under similarity transformations,
algebraic identities, and exact agreement between independent algorithms
(MPC without active constraints equals LQR, MHE without active
constraints equals the Kalman filter, H-infinity norms equal dense
frequency sweeps). A separate GitHub Actions workflow compares Riccati
and Lyapunov solutions, norms, Hankel singular values, discretization,
pole placement and H-infinity synthesis against python-control and
SLICOT.

## Related packages

controlcompR belongs to a family of R packages by the same author:
[symcompR](https://github.com/adamrsophie/symcompR) (symbolic
computation), [modelcompR](https://github.com/adamrsophie/modelcompR)
(equation-based modeling), [jaxR](https://github.com/adamrsophie/jaxR)
(differentiable computation), [z3r](https://github.com/adamrsophie/z3r)
(constraint solving),
[optcompR](https://github.com/adamrsophie/optcompR) (optimization),
[graphlearnR](https://github.com/adamrsophie/graphlearnR) (graph
learning) and [pdecompR](https://github.com/adamrsophie/pdecompR)
(partial differential equations). controlcompR works on its own.

## License

GPL-3.
