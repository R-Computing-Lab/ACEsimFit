test_that("Power and sample size can match1", {
  expect_equal(round(Power_LS(N1 = 57, N2 = 114, h2 = .5, c2 = .2, R1 = 1, R2 = .50),1),.8)
})

test_that("Power and sample size can match2", {
     expect_equal(round(Power_LS(p_N1 = .333333, power = .8,h2 = .5, c2 = .2, R1 = 1, R2 = .50)[1]),57)
})
