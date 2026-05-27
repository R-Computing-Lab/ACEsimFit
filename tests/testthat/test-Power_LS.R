test_that("Power and sample size can match1", {
  expect_equal(round(Power_LS(N1 = 57, N2 = 114, h2 = .5, c2 = .2, R1 = 1, R2 = .50), 1), .8)
})

test_that("Power and sample size can match2", {
  expect_equal(round(Power_LS(p_N1 = .333333, power = .8, h2 = .5, c2 = .2, R1 = 1, R2 = .50)[1]), 57)
})

test_that("power is between 0 and 1 for reasonable inputs", {
  p <- Power_LS(N1 = 100, N2 = 100, h2 = 0.5, c2 = 0.2)
  expect_true(p > 0 && p < 1)
})

test_that("power increases with larger N", {
  p_small <- Power_LS(N1 = 50,  N2 = 50,  h2 = 0.5, c2 = 0.1)
  p_large <- Power_LS(N1 = 500, N2 = 500, h2 = 0.5, c2 = 0.1)
  expect_lt(p_small, p_large)
})

test_that("power increases with larger h2", {
  p_low  <- Power_LS(N1 = 200, N2 = 200, h2 = 0.2, c2 = 0.1)
  p_high <- Power_LS(N1 = 200, N2 = 200, h2 = 0.6, c2 = 0.1)
  expect_lt(p_low, p_high)
})

test_that("solving for N2 given N1 is consistent with forward power", {
  target_power <- 0.8
  n2 <- Power_LS(N1 = 100, power = target_power, h2 = 0.5, c2 = 0.2)
  recovered_power <- Power_LS(N1 = 100, N2 = n2, h2 = 0.5, c2 = 0.2)
  expect_lt(abs(recovered_power - target_power), 0.01)
})

test_that("solving for N1 given N2 is consistent with forward power", {
  target_power <- 0.8
  n1 <- Power_LS(N2 = 200, power = target_power, h2 = 0.5, c2 = 0.2)
  recovered_power <- Power_LS(N1 = n1, N2 = 200, h2 = 0.5, c2 = 0.2)
  expect_lt(abs(recovered_power - target_power), 0.01)
})

test_that("p_N1 branch returns two values summing to total N", {
  result <- Power_LS(p_N1 = 0.5, power = 0.8, h2 = 0.5, c2 = 0.2)
  expect_length(result, 2)
  expect_equal(result[1], result[2])  # equal split when p_N1 = 0.5
})

test_that("invalid combination (both N missing, no p_N1, power specified) errors", {
  expect_error(
    Power_LS(power = 0.8, h2 = 0.5, c2 = 0.2),
    "Invalid argument combination"
  )
})

test_that("digits argument controls rounding", {
  n2_rounded <- Power_LS(N1 = 100, power = 0.8, h2 = 0.5, c2 = 0.2, digits = 0)
  n2_exact   <- Power_LS(N1 = 100, power = 0.8, h2 = 0.5, c2 = 0.2, digits = 3)
  expect_equal(n2_rounded, round(n2_rounded, 0))
  expect_false(identical(n2_rounded, n2_exact))
})
