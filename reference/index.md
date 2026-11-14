# Package index

## Models

- [`ctrl_ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ss.md)
  : State-space models
- [`ctrl_tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf.md)
  : Transfer-function models
- [`ctrl_zpk()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_zpk.md)
  : Zero-pole-gain models
- [`ctrl_gain()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_gain.md)
  : Static gain model
- [`ctrl_tf2ss()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
  [`ctrl_ss2tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
  [`ctrl_zpk2tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
  [`ctrl_tf2zpk()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf2ss.md)
  : Conversions between model representations
- [`ctrl_ssdata()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ssdata.md)
  [`ctrl_dims()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ssdata.md)
  [`ctrl_is_discrete()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ssdata.md)
  [`ctrl_is_continuous()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ssdata.md)
  [`ctrl_sample_time()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ssdata.md)
  : State-space data and dimensions
- [`ctrl_transform()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_transform.md)
  : Similarity transformation of a state-space model
- [`ctrl_signal_names()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_signal_names.md)
  : Signal names of a model
- [`ctrl_is_proper()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_is_proper.md)
  : Test whether a transfer function is proper
- [`is_ctrl_ss()`](https://adamrsophie.github.io/controlcompR/reference/is_ctrl_ss.md)
  [`is_ctrl_tf()`](https://adamrsophie.github.io/controlcompR/reference/is_ctrl_ss.md)
  [`is_ctrl_system()`](https://adamrsophie.github.io/controlcompR/reference/is_ctrl_ss.md)
  : Test for control system objects
- [`` `[`( ``*`<ctrl_ss>`*`)`](https://adamrsophie.github.io/controlcompR/reference/sub-.ctrl_ss.md)
  : Select inputs and outputs of a model
- [`summary(`*`<ctrl_ss>`*`)`](https://adamrsophie.github.io/controlcompR/reference/summary.ctrl_ss.md)
  : Summaries of control models
- [`ctrl_c2d()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_c2d.md)
  : Discretize a continuous time model
- [`ctrl_d2c()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_d2c.md)
  : Convert a discrete time model to continuous time

## Interconnections

- [`ctrl_series()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_series.md)
  [`ctrl_parallel()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_series.md)
  [`ctrl_feedback()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_series.md)
  [`ctrl_append()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_series.md)
  : Series, parallel and feedback connections
- [`Ops(`*`<ctrl_system>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_arith.md)
  : Arithmetic on models
- [`ctrl_inv()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_inv.md)
  : Inverse of a model
- [`ctrl_connect()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_connect.md)
  : Connect models by signal names
- [`ctrl_sumblock()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_sumblock.md)
  : Summing junction
- [`ctrl_closed_loop()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_closed_loop.md)
  : Closed-loop systems
- [`ctrl_sensitivity()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_sensitivity.md)
  [`ctrl_gangof4()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_sensitivity.md)
  : Sensitivity functions
- [`ctrl_lft()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lft.md)
  [`ctrl_plant_partition()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lft.md)
  : Generalized plants and linear fractional transformations

## Structure and stability

- [`ctrl_poles()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_poles.md)
  : Poles of a model
- [`ctrl_zeros()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_zeros.md)
  : Invariant zeros
- [`ctrl_damp()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_damp.md)
  : Natural frequencies and damping ratios
- [`ctrl_dcgain()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_dcgain.md)
  : Steady-state (DC) gain
- [`ctrl_is_stable()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_is_stable.md)
  : Test asymptotic stability
- [`ctrl_stability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_stability.md)
  : Stability analysis
- [`ctrl_controllability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_controllability.md)
  : Controllability and stabilizability analysis
- [`ctrl_observability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_observability.md)
  : Observability and detectability analysis
- [`ctrl_ctrb()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ctrb.md)
  [`ctrl_obsv()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ctrb.md)
  : Controllability and observability matrices
- [`ctrl_is_stabilizable()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_is_stabilizable.md)
  [`ctrl_is_detectable()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_is_stabilizable.md)
  [`ctrl_is_controllable()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_is_stabilizable.md)
  [`ctrl_is_observable()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_is_stabilizable.md)
  : Stabilizability and detectability tests
- [`ctrl_gram()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_gram.md)
  : Controllability and observability Gramians
- [`ctrl_minreal()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_minreal.md)
  : Minimal realization
- [`ctrl_tf_minreal()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_tf_minreal.md)
  : Pole-zero cancellation for transfer functions
- [`ctrl_stabsep()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_stabsep.md)
  : Stable and unstable decomposition

## Matrix equations

- [`ctrl_lyap()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lyap.md)
  : Lyapunov and Stein equations
- [`ctrl_sylvester()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_sylvester.md)
  : Sylvester equation
- [`ctrl_care()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_care.md)
  [`ctrl_dare()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_care.md)
  : Algebraic Riccati equations

## Time responses

- [`ctrl_simulate()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_simulate.md)
  : Simulate a linear model
- [`ctrl_step()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md)
  [`ctrl_impulse()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md)
  [`ctrl_initial()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step.md)
  : Step, impulse and initial-condition responses
- [`ctrl_step_info()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_step_info.md)
  : Step-response characteristics
- [`as.data.frame(`*`<ctrl_response>`*`)`](https://adamrsophie.github.io/controlcompR/reference/as.data.frame.ctrl_response.md)
  : Long data frame of a time response

## Frequency responses, margins and norms

- [`ctrl_freqresp()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_freqresp.md)
  : Frequency response
- [`ctrl_freq_grid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_freq_grid.md)
  : Frequency grid for a model
- [`ctrl_bode()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_bode.md)
  [`ctrl_nyquist()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_bode.md)
  [`ctrl_nichols()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_bode.md)
  : Bode, Nyquist and Nichols data
- [`ctrl_margin()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_margin.md)
  : Gain, phase and delay margins
- [`ctrl_disk_margin()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_disk_margin.md)
  : Disk margins
- [`ctrl_bandwidth()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_bandwidth.md)
  [`ctrl_peak()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_bandwidth.md)
  : Bandwidth and resonant peak
- [`ctrl_norm()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_norm.md)
  : System norms
- [`ctrl_rlocus()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_rlocus.md)
  : Root locus
- [`ctrl_pzmap()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pzmap.md)
  : Pole-zero map data
- [`plot(`*`<ctrl_kalman_filter>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_idcompare>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_resid>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_mpc_simulation>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_nl_response>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_pid_simulation>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_response>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_bode>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_nyquist>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_nichols>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_pzmap>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  [`plot(`*`<ctrl_rlocus>`*`)`](https://adamrsophie.github.io/controlcompR/reference/ctrl_plots.md)
  : Plot methods

## State feedback and estimation

- [`ctrl_lqr()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lqr.md)
  : Linear quadratic regulator
- [`ctrl_dlqr_finite()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_dlqr_finite.md)
  : Finite-horizon discrete LQR
- [`ctrl_reference_gain()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_reference_gain.md)
  : Reference (feedforward) gain for state feedback
- [`ctrl_place()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_place.md)
  : Pole placement
- [`ctrl_observer()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_observer.md)
  : Luenberger observer
- [`ctrl_kalman()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman.md)
  : Steady-state Kalman filter
- [`ctrl_kalman_filter()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman_filter.md)
  [`ctrl_kalman_smoother()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_kalman_filter.md)
  : Time-varying discrete Kalman filter and smoother
- [`ctrl_lqg()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_lqg.md)
  : LQG controller
- [`ctrl_simulate_lqg()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_simulate_lqg.md)
  : Stochastic closed-loop simulation of an LQG controller

## PID control

- [`ctrl_pid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid.md)
  : PID controllers
- [`ctrl_pid_tf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_tf.md)
  : Linear transfer functions of a PID controller
- [`ctrl_pid_update()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_update.md)
  : One step of a discrete PID controller
- [`ctrl_simulate_pid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_simulate_pid.md)
  : Closed-loop simulation with a discrete PID controller
- [`ctrl_pid_tune()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_pid_tune.md)
  : PID tuning rules
- [`ctrl_fit_fopdt()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_fit_fopdt.md)
  : Fit a first-order plus dead-time model to a step response

## Performance assessment

- [`ctrl_performance()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_performance.md)
  : Tracking and control performance metrics
- [`ctrl_compare_controllers()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_compare_controllers.md)
  : Compare controllers on a common scenario

## Model predictive control and moving horizon estimation

- [`ctrl_mpc()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mpc.md)
  : Linear model predictive control
- [`ctrl_mpc_step()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mpc_step.md)
  : One MPC step
- [`ctrl_mpc_qp()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mpc_qp.md)
  : Inspect the MPC quadratic program
- [`ctrl_mpc_estimator()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mpc_estimator.md)
  : Offset-free MPC estimator
- [`ctrl_simulate_mpc()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_simulate_mpc.md)
  : Closed-loop MPC simulation
- [`ctrl_mhe()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mhe.md)
  [`ctrl_mhe_run()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_mhe.md)
  : Moving horizon estimation for linear models

## Robust control

- [`ctrl_h2syn()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_h2syn.md)
  : H2 optimal output-feedback synthesis
- [`ctrl_hinfsyn()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_hinfsyn.md)
  : H-infinity output-feedback synthesis (experimental)
- [`ctrl_augw()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_augw.md)
  : Mixed-sensitivity generalized plant
- [`ctrl_robust_stability()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_robust_stability.md)
  : Robust stability under unstructured uncertainty
- [`ctrl_robust_performance()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_robust_performance.md)
  : Robust performance under multiplicative output uncertainty

## Model reduction

- [`ctrl_hsvd()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_hsvd.md)
  : Hankel singular values
- [`ctrl_balreal()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_balreal.md)
  [`ctrl_balred()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_balreal.md)
  : Balanced realization and balanced truncation

## System identification

- [`ctrl_iddata()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_iddata.md)
  : Identification data

- [`ctrl_idsplit()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_idsplit.md)
  : Split identification data into estimation and validation parts

- [`ctrl_arx()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_arx.md)
  : ARX model estimation

- [`ctrl_arx_order()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_arx_order.md)
  : ARX model order selection

- [`ctrl_n4sid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_n4sid.md)
  : Subspace state-space identification (PO-MOESP)

- [`ctrl_from_sysid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_from_sysid.md)
  :

  Convert a `sysid` model

- [`ctrl_compare()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_compare.md)
  : Compare identified models with data

- [`ctrl_resid()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_resid.md)
  : Residual analysis for identified models

## Nonlinear systems

- [`ctrl_nlsys()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nlsys.md)
  : Nonlinear state-space models
- [`ctrl_nl_simulate()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nl_simulate.md)
  : Simulate a nonlinear model
- [`ctrl_nl_discretize()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nl_discretize.md)
  : Discretize a continuous nonlinear model
- [`ctrl_linearize()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_linearize.md)
  : Linearize a nonlinear model
- [`ctrl_equilibrium()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_equilibrium.md)
  : Equilibrium points
- [`ctrl_ekf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ekf.md)
  [`ctrl_ukf()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_ekf.md)
  : Extended and unscented Kalman filters
- [`ctrl_nmpc()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nmpc.md)
  [`ctrl_nmpc_step()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nmpc.md)
  [`ctrl_simulate_nmpc()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nmpc.md)
  : Nonlinear model predictive control (experimental)
- [`ctrl_nl_trajopt()`](https://adamrsophie.github.io/controlcompR/reference/ctrl_nl_trajopt.md)
  : Nonlinear trajectory optimization
