test_that("Results Level Check", {
  expect_equal(length(Sim_Fit(nIter = 5, saveRaw = TRUE)), 5)
})

test_that("Results Category Check", {
  expect_equal(length(Sim_Fit(nIter = 5, saveRaw = TRUE)[[1]]), 2)
})

test_that("each iteration has Results and data slots", {
  out <- Sim_Fit(nIter = 2, saveRaw = TRUE)
  expect_named(out[[1]], c("Results", "data"))
})

test_that("saveRaw = FALSE stores NA in data slot", {
  out <- Sim_Fit(nIter = 2, saveRaw = FALSE)
  expect_true(is.na(out[[1]]$data))
})

test_that("saveRaw = TRUE stores a data.frame in data slot", {
  out <- Sim_Fit(nIter = 2, saveRaw = TRUE)
  expect_s3_class(out[[1]]$data, "data.frame")
})

test_that("iteration names follow 'Iteration<i>' pattern", {
  out <- Sim_Fit(nIter = 3, saveRaw = FALSE)
  expect_equal(names(out), c("Iteration1", "Iteration2", "Iteration3"))
})

test_that("Results slot contains nest and summary elements", {
  out <- Sim_Fit(nIter = 2, saveRaw = FALSE)
  expect_true(all(c("nest", "summary") %in% names(out[[1]]$Results)))
})

test_that("SSeed makes results reproducible", {
  out1 <- Sim_Fit(nIter = 2, SSeed = 42, saveRaw = TRUE)
  out2 <- Sim_Fit(nIter = 2, SSeed = 42, saveRaw = TRUE)
  expect_equal(out1[[1]]$data, out2[[1]]$data)
})

test_that("different SSeed produces different data", {
  out1 <- Sim_Fit(nIter = 2, SSeed = 1,  saveRaw = TRUE)
  out2 <- Sim_Fit(nIter = 2, SSeed = 99, saveRaw = TRUE)
  expect_false(identical(out1[[1]]$data, out2[[1]]$data))
})
