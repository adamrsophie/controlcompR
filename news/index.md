# Changelog

## controlcompR 0.1.0

- First version.
  - Models: continuous and discrete state-space, transfer-function,
    zero-pole-gain and static gain models with MIMO support and signal
    names; conversions; zero-order hold, first-order hold, Tustin and
    matched discretization and
    [`ctrl_d2c()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_d2c.md).
  - Interconnections: series, parallel, feedback with well-posedness
    checks, block diagrams by signal name, closed loops, sensitivity
    functions and linear fractional transformations.
  - Analysis: poles, transmission zeros, damping, stability with
    explicit tolerances, staircase controllability and observability,
    Gramians, minimal realizations, stable and unstable decomposition.
  - Matrix equations: Lyapunov, Stein and Sylvester solvers; continuous
    and discrete Riccati solvers based on ordered Schur forms, with
    residual and stability verification.
  - Time and frequency responses, step metrics, gain, phase, delay and
    disk margins, bandwidth, H2 and H-infinity norms, root locus.
  - Design: LQR (continuous, discrete, finite horizon), robust pole
    placement, observers, Kalman filters and the RTS smoother, LQG, PID
    controllers with filters, anti-windup and bumpless transfer, tuning
    rules.
  - Linear MPC with OSQP: input, rate, state, output and terminal
    constraints, soft constraints, offset-free tracking and explicit
    failure handling; moving horizon estimation with warm starts.
  - Robust control: H2 synthesis and experimental, certified H-infinity
    synthesis for continuous time with `D11 = 0`; mixed-sensitivity
    plants; robust stability and robust performance tests for
    unstructured uncertainty.
  - Model reduction by balanced truncation and residualization.
  - System identification: ARX, order selection, subspace
    identification, validation tools and conversion of `sysid` models.
  - Nonlinear systems: simulation, linearization, equilibria,
    discretization, EKF, UKF, experimental nonlinear MPC with adjoint
    gradients and trajectory optimization with terminal constraints.
  - Fifteen articles and a reference comparison against python-control
    and SLICOT.
- Numerical contracts worth knowing:
  - Balanced truncation checks its error bound on a frequency sample and
    stops with an error when the bound cannot be certified.
  - [`ctrl_robust_performance()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_robust_performance.md)
    samples its peaks (including zero and infinite frequency) and
    reports a separate `certified` flag from a rigorous upper bound; the
    SISO exactness flag is `criterion_exact`.
  - Nyquist encirclement counts are `NA` when the frequency sweep cannot
    resolve the contour.
  - [`ctrl_zeros()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_zeros.md)
    stops with an error for systems whose pencil is singular (no full
    normal rank).
  - Step, impulse and initial responses are applied at time zero; grids
    that start later report the same response at their own instants.
  - [`ctrl_step_info()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step_info.md)
    reports `final_value_known` separately from `settled`, which is now
    always based on the simulated horizon.
  - [`ctrl_pid_update()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_update.md)
    has `init = "hold"` (default) and `init = "zero"` (zero history,
    matching the zero-state response of
    [`ctrl_pid_tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_tf.md)).
  - EKF and UKF results have class `ctrl_nl_filter` and are not accepted
    by
    [`ctrl_kalman_smoother()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman_filter.md);
    infinite measurements are rejected (only `NA` marks missing data).
  - NMPC simulation needs a discrete plant with the controller’s sample
    interval; the `"hold"` fallback is clipped to the input bounds and
    the first held input is set by `u0`.
  - Numeric scalars scale every channel of MIMO models.
