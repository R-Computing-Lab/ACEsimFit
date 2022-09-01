test_that("two group exist", {
  expect_equal(length(unique(kinsim_double(GroupRel = c(.98,.644),
                                           ace2 = c(5,7,9),
                                           ifComb = TRUE)$GroupName)), 2)
})

test_that("nrow correct", {
     expect_equal(nrow(kinsim_double(GroupSizes = c(59,131),
                                     GroupRel = c(1,.644),
                                     ace1 = c(.5,.7,.19),
                                     ifComb = TRUE)), 190)
})

test_that("ncol correct", {
     expect_equal(ncol(kinsim_double(GroupSizes = c(159,231),
                                     GroupRel = c(.98,.5),
                                     ace2 = c(5,7,9),
                                     ifComb = FALSE)), 12)
})
