
test_that("step info handles direct feedthrough", {
  si <- ctrl_step_info(ctrl_tf(c(1, 2), c(1, 1)))
  expect_equal(si$final_value, 2)
  expect_false(is.na(si$settling_time))
  Gd <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.5, method = "tustin")
  expect_equal(ctrl_step_info(Gd)$final_value, 1, tolerance = 1e-10)
})

test_that("performance metrics check shapes and per-input limits", {
  t <- 0:1
  p <- ctrl_performance(t, c(0, 0), 0, u = rbind(c(0.5, 5), c(0.5, 5)),
                        u_min = -c(1, 10), u_max = c(1, 10))
  expect_equal(p$max_violation, 0)
  expect_error(ctrl_performance(0:5, 1:3), "rows")
  plant <- ctrl_c2d(ctrl_tf(1, c(1, 1)), 0.1)
  mpc <- ctrl_mpc(plant, 5, Q = 1, R = 0.1, u_min = -1, u_max = 1)
  sim <- ctrl_simulate_mpc(mpc, 1, 20)
  expect_silent(ctrl_performance(sim$time, sim$y, 0, u = sim$u))
})
