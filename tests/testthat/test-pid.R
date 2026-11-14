test_that("PID transfer function matches the parallel form", {
  pid <- ctrl_pid(kp = 2, ki = 3, kd = 0.5, Tf = 0.05)
  C <- ctrl_pid_tf(pid)
  w <- 1.7
  s <- 1i * w
  want <- 2 + 3 / s + 0.5 * s / (0.05 * s + 1)
  got <- ctrl_freqresp(C, w)$response[1, 1, 1]
  expect_equal(got, want, tolerance = 1e-10)
  expect_error(ctrl_pid(kp = 1, kd = 1), "filter")
})

test_that("discrete PID update matches its transfer function when unsaturated", {
  pid <- ctrl_pid(kp = 1.5, ki = 0.8, kd = 0.3, Tf = 0.1, dt = 0.05, b = 1, c = 1)
  set.seed(2)
  e <- rnorm(40)
  st <- NULL
  u <- numeric(40)
  for (k in 1:40) {
    out <- ctrl_pid_update(pid, st, r = e[k], y = 0)
    st <- out$state
    u[k] <- out$u
  }
  C <- ctrl_pid_tf(pid)
  ref <- ctrl_simulate(C, input = e, times = (0:39) * 0.05)$y[, 1, 1]
  # the update starts the derivative from the first sample (no kick), so it
  # differs from the zero-state transfer function by a decaying term
  a <- pid$Tf / (pid$Tf + 0.05)
  g <- pid$kd / (pid$Tf + 0.05)
  expect_equal(ref - u, g * e[1] * a^(0:39), tolerance = 1e-8)
})

test_that("anti-windup limits integrator growth under saturation", {
  plant <- ctrl_tf(1, c(1, 1))
  mk <- function(aw) ctrl_pid(kp = 1, ki = 2, u_min = -0.5, u_max = 0.5,
                              antiwindup = aw, dt = 0.05)
  sim_none <- ctrl_simulate_pid(plant, mk("none"), reference = 1, steps = 200)
  sim_bc <- ctrl_simulate_pid(plant, mk("backcalc"), reference = 1, steps = 200)
  sim_cl <- ctrl_simulate_pid(plant, mk("clamp"), reference = 1, steps = 200)
  expect_true(all(abs(sim_bc$u) <= 0.5 + 1e-12))
  expect_gt(max(sim_none$I), max(sim_bc$I))
  expect_gt(max(sim_none$I), max(sim_cl$I))
})

test_that("bumpless transfer keeps the output continuous", {
  plant <- ctrl_tf(1, c(1, 1))
  pid <- ctrl_pid(kp = 1, ki = 1, dt = 0.1)
  manual <- c(rep(0.3, 30), rep(NA, 30))
  sim <- ctrl_simulate_pid(plant, pid, reference = 0.3, steps = 60, manual = manual)
  expect_lt(abs(sim$u[31] - sim$u[30]), 0.05)
})

test_that("tuning rules return sensible controllers", {
  plant <- ctrl_tf(1, c(1, 3, 3, 1))
  zn <- ctrl_pid_tune("ziegler_nichols", plant = plant, derivative_filter = 10)
  expect_equal(zn$kp, 0.6 * 8, tolerance = 1e-6)
  fit <- ctrl_fit_fopdt(plant)
  expect_equal(unname(fit["K"]), 1, tolerance = 1e-3)
  simc <- ctrl_pid_tune("simc", fopdt = fit)
  expect_gt(simc$kp, 0)
  amigo <- ctrl_pid_tune("amigo", fopdt = fit)
  expect_gt(amigo$kp, 0)
})

test_that("step info and performance metrics follow their definitions", {
  G <- ctrl_tf(1, c(1, 1))
  info <- ctrl_step_info(G)
  expect_equal(info$rise_time, log(9), tolerance = 1e-2)
  expect_equal(info$settling_time, -log(0.02), tolerance = 1e-2)
  expect_equal(info$overshoot, 0)
  w <- 2; z <- 0.3
  G2 <- ctrl_tf(w^2, c(1, 2 * z * w, w^2))
  info2 <- ctrl_step_info(G2)
  expect_equal(info2$overshoot, 100 * exp(-pi * z / sqrt(1 - z^2)), tolerance = 1e-2)
  t <- seq(0, 20, by = 0.001)
  perf <- ctrl_performance(t, 1 - exp(-t), 1)
  expect_equal(perf$IAE, 1, tolerance = 1e-4)
  expect_equal(perf$ISE, 0.5, tolerance = 1e-4)
})

test_that("manual to automatic transfer is bumpless for every integral method", {
  for (im in c("forward", "backward", "tustin")) {
    pid <- ctrl_pid(kp = 1, ki = 1, kd = 0.5, derivative_filter = 5, dt = 1,
                    integral_method = im)
    o1 <- ctrl_pid_update(pid, NULL, r = 1, y = 0, manual = 2)
    o2 <- ctrl_pid_update(pid, o1$state, r = 1, y = 0)
    expect_equal(o2$u, 2, tolerance = 1e-12)
  }
})

test_that("PID simulation refuses plants with feedthrough", {
  expect_error(ctrl_simulate_pid(ctrl_gain(1, dt = 1), ctrl_pid(kp = 1, dt = 1), steps = 3),
               "feedthrough")
})
