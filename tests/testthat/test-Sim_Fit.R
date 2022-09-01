test_that("Results Level Check", {
  expect_equal(length(Sim_Fit(nIter = 5,saveRaw=TRUE)), 5)
})

test_that("Results Category Check", {
     expect_equal(length(Sim_Fit(nIter = 5,saveRaw=TRUE)[[1]]), 2)
})

