test_that("single simulation works", {
     expect_equal(ncol(kinsim_single(
          name = "testtesttest",
          Rel = .8,
          r_c = .98,
          n = 1000,
          mu = 2,
          ace = c(2,2,6)
     )), 12)
})
