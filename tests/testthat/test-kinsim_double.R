test_that("two group exist", {
  expect_equal(length(unique(kinsim_double(
    GroupRel = c(.98, .644),
    ace2 = c(5, 7, 9),
    ifComb = TRUE
  )$GroupName)), 2)
})

test_that("nrow correct", {
  expect_equal(nrow(kinsim_double(
    GroupSizes = c(59, 131),
    GroupRel = c(1, .644),
    ace1 = c(.5, .7, .19),
    ifComb = TRUE
  )), 190)
})

test_that("ncol correct", {
  expect_equal(ncol(kinsim_double(
    GroupSizes = c(159, 231),
    GroupRel = c(.98, .5),
    ace2 = c(5, 7, 9),
    ifComb = FALSE
  )), 12)
})

test_that("direct approach: group sizes match GroupSizes", {
  df <- kinsim_double(GroupSizes = c(80, 120), ifComb = FALSE)
  counts <- table(df$GroupName)
  expect_equal(as.integer(counts[["KinPair1"]]), 80)
  expect_equal(as.integer(counts[["KinPair2"]]), 120)
})

test_that("column names are correct", {
  df <- kinsim_double()
  expect_named(df, c("GroupName", "R", "r_c", "id", "A1", "A2", "C1", "C2", "E1", "E2", "y1", "y2"))
})

test_that("R column values match GroupRel in direct approach", {
  df <- kinsim_double(GroupRel = c(1, 0.5), ifComb = FALSE)
  r_vals <- sort(unique(df$R))
  expect_equal(r_vals, c(0.5, 1.0))
})

test_that("combination approach: both extremes (both Rel in {1, .5})", {
  df <- kinsim_double(
    GroupRel = c(1, 0.5),
    GroupSizes = c(100, 100),
    ifComb = TRUE
  )
  expect_equal(nrow(df), 200)
  expect_equal(length(unique(df$GroupName)), 2)
})

test_that("combination approach: both Rel between .5 and 1", {
  df <- kinsim_double(
    GroupRel = c(0.75, 0.65),
    GroupSizes = c(100, 100),
    ifComb = TRUE
  )
  expect_equal(nrow(df), 200)
  # R column should be set to the target Rel value for both groups
  expect_true(all(df$R[df$GroupName == "KinPair1"] == 0.75))
  expect_true(all(df$R[df$GroupName == "KinPair2"] == 0.65))
})

test_that("combination approach: group1 Rel in {1,.5}, group2 Rel is intermediate", {
  df <- kinsim_double(
    GroupRel = c(1, 0.75),
    GroupSizes = c(100, 100),
    ifComb = TRUE
  )
  expect_equal(nrow(df), 200)
  expect_true(all(df$R[df$GroupName == "KinPair2"] == 0.75))
})

test_that("combination approach: group1 Rel is intermediate, group2 Rel in {1,.5}", {
  df <- kinsim_double(
    GroupRel = c(0.75, 0.5),
    GroupSizes = c(100, 100),
    ifComb = TRUE
  )
  expect_equal(nrow(df), 200)
  expect_true(all(df$R[df$GroupName == "KinPair1"] == 0.75))
})
