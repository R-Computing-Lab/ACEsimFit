test_that("single simulation works", {
  expect_equal(ncol(kinsim_single(
    name = "testtesttest",
    Rel = .8,
    r_c = .98,
    n = 1000,
    mu = 2,
    ace = c(2, 2, 6)
  )), 12)
})

test_that("output has correct number of rows", {
  df <- kinsim_single(n = 150)
  expect_equal(nrow(df), 150)
})

test_that("column names are correct", {
  df <- kinsim_single()
  expect_named(df, c("GroupName", "R", "r_c", "id", "A1", "A2", "C1", "C2", "E1", "E2", "y1", "y2"))
})

test_that("GroupName column matches name argument", {
  df <- kinsim_single(name = "Twins_MZ", n = 50)
  expect_true(all(df$GroupName == "Twins_MZ"))
})

test_that("R column matches Rel argument", {
  df <- kinsim_single(Rel = 0.5, n = 50)
  expect_true(all(df$R == 0.5))
})

test_that("r_c column matches r_c argument", {
  df <- kinsim_single(r_c = 0.75, n = 50)
  expect_true(all(df$r_c == 0.75))
})

test_that("id column is 1:n without gaps", {
  df <- kinsim_single(n = 80)
  expect_equal(df$id, seq_len(80))
})

test_that("y values are approximately centred on mu", {
  set.seed(42)
  df <- kinsim_single(n = 5000, mu = 10, ace = c(0.5, 0.3, 0.2))
  expect_lt(abs(mean(c(df$y1, df$y2)) - 10), 0.2)
})

test_that("total variance is approximately sum of ACE components", {
  set.seed(7)
  df <- kinsim_single(n = 10000, mu = 0, ace = c(0.4, 0.3, 0.3))
  observed_var <- var(c(df$y1, df$y2))
  expect_lt(abs(observed_var - 1.0), 0.1)
})

test_that("MZ pair correlation is approximately A + C for r_c = 1, Rel = 1", {
  set.seed(99)
  df <- kinsim_single(n = 5000, Rel = 1, r_c = 1, ace = c(0.6, 0.2, 0.2))
  pair_cor <- cor(df$y1, df$y2)
  # Expected: (A + C) / (A + C + E) = 0.8
  expect_lt(abs(pair_cor - 0.8), 0.05)
})

test_that("DZ pair correlation is approximately 0.5*A + C for r_c = 1, Rel = 0.5", {
  set.seed(123)
  df <- kinsim_single(n = 5000, Rel = 0.5, r_c = 1, ace = c(0.6, 0.2, 0.2))
  pair_cor <- cor(df$y1, df$y2)
  # Expected: (0.5*0.6 + 0.2) / 1 = 0.5
  expect_lt(abs(pair_cor - 0.5), 0.05)
})
