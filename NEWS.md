# controlcompR 0.1.0

* First version.
  * Models: continuous and discrete state-space, transfer-function,
    zero-pole-gain and static gain models with MIMO support and signal
    names; conversions; zero-order hold, first-order hold, Tustin and
    matched discretization and `ctrl_d2c()`.
  * Interconnections: series, parallel, feedback with well-posedness
    checks, block diagrams by signal name, closed loops, sensitivity
    functions and linear fractional transformations.
  * Analysis: poles, transmission zeros, damping, stability with explicit
    tolerances, staircase controllability and observability, Gramians,
    minimal realizations, stable and unstable decomposition.
  * Matrix equations: Lyapunov, Stein and Sylvester solvers; continuous
    and discrete Riccati solvers based on ordered Schur forms, with
    residual and stability verification.
  * Time and frequency responses, step metrics, gain, phase, delay and
    disk margins, bandwidth, H2 and H-infinity norms, root locus.
  * Design: LQR (continuous, discrete, finite horizon), robust pole
    placement, observers, Kalman filters and the RTS smoother, LQG, PID
    controllers with filters, anti-windup and bumpless transfer, tuning
    rules.
  * Linear MPC with OSQP: input, rate, state, output and terminal
    constraints, soft constraints, offset-free tracking and explicit
    failure handling; moving horizon estimation with warm starts.
  * Robust control: H2 synthesis and experimental, certified H-infinity
    synthesis for continuous time with `D11 = 0`; mixed-sensitivity
    plants; robust stability and robust performance tests for unstructured
    uncertainty.
  * Model reduction by balanced truncation and residualization.
  * System identification: ARX, order selection, subspace identification,
    validation tools and conversion of `sysid` models.
  * Nonlinear systems: simulation, linearization, equilibria,
    discretization, EKF, UKF, experimental nonlinear MPC with adjoint
    gradients and trajectory optimization with terminal constraints.
  * Fifteen articles and a reference comparison against python-control
    and SLICOT.
* Numerical contracts worth knowing:
  * Gramian factors for Hankel singular values, balanced realization,
    balanced truncation and the H2 norm are computed directly with
    Hammarling's method, so the Gramians are never formed and rotated or
    badly scaled realizations keep their weak directions.
  * Balanced truncation checks its error bound on a frequency sample and
    against the error system's H-infinity norm, stops with an error when
    its samples contradict the bound, and reports `bound_status`
    (`"certified"` or `"estimated"`) and `error_norm`. `ctrl_balreal()`
    runs the same check on the removed states.
  * Stability is decided on a realization balanced together with `B` and
    `C`, without a fixed unit floor in continuous time, and with per-pole
    bands from the first-order eigenvalue error bound (`pole_tol`).
  * Covariances are inverted on their support after scaling to unit
    diagonal in the UKF update, the RTS smoother and the MHE arrival
    cost, and every direction of a covariance accepted as positive
    definite is kept, also when rotated away from the axes; the MHE QP is
    solved in the prior's own coordinates scaled to unit curvature.
  * Loop internal stability checks use the same policy as
    `ctrl_stability()`.
  * Robust stability decisions and H-infinity synthesis certification use
    the upper bound of the norm computation (`upper_bound`); an attained
    zero H-infinity optimum is returned with `gamma = 0`.
  * `ctrl_robust_performance()` samples its peaks (including zero and
    infinite frequency) and reports a separate `certified` flag from a
    rigorous upper bound; the SISO exactness flag is `criterion_exact`.
  * Nyquist encirclement counts come from the argument principle (unstable
    closed-loop minus open-loop roots) and are `NA` when a closed-loop pole
    lies on the stability boundary or `1 + D = 0`.
  * `ctrl_zeros()` stops with an error for systems whose pencil is
    singular (no full normal rank), decided by evaluating the transfer
    matrix away from every computed eigenvalue and pole. With a
    well-conditioned `D`, zeros too large for the pencil to separate from
    infinity are taken from `A - B D^-1 C`; otherwise an eigenvalue that
    cannot be told apart from infinity gives a warning.
  * Step, impulse and initial responses are applied at time zero; grids
    that start later report the same response at their own instants.
  * `ctrl_step_info()` reports `final_value_known` separately from
    `settled`, which is now always based on the simulated horizon.
  * `ctrl_pid_update()` has `init = "hold"` (default) and `init = "zero"`
    (zero history, matching the zero-state response of `ctrl_pid_tf()`).
  * EKF and UKF results have class `ctrl_nl_filter` and are not accepted
    by `ctrl_kalman_smoother()`; infinite measurements are rejected (only
    `NA` marks missing data).
  * NMPC simulation needs a discrete plant with the controller's sample
    interval; the `"hold"` fallback is clipped to the input bounds and
    the first held input is set by `u0`.
  * Numeric scalars scale every channel of MIMO models; dividing a model
    by a square numeric matrix multiplies by its inverse on the input side.
  * Riccati solvers equilibrate the equation (weight scale, state units,
    time scale and input units) before solving, so solutions are
    invariant to the units of states, inputs, costs and time; the
    reported residual is that of the equilibrated equation.
  * The Kalman filter, EKF, UKF and MHE reject non-numeric data and rows
    that are only partly missing; an entirely `NA` row is a missing
    sample.
  * Fast or slow SISO transfer functions are realized with scaled states
    when the scaled coefficients stay representable; otherwise the
    unscaled companion form is used.
  * Matched discretization treats a root as an origin root only when it
    is negligible both against the other roots and against the sample
    interval.
  * The PID loop through plant feedthrough uses the exact affine slope of
    the controller output.
  * Scalar model arithmetic that overflows, nonfinite Kalman prior means
    and negative NMPC penalties are rejected.
  * Static loops adopt the clock of a discrete weight in
    `ctrl_robust_performance()`; synthesis rank checks on `D12` and `D21`
    are relative to each block.
