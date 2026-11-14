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
