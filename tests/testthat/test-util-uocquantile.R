test_that("uocquantile ordered-factor branch matches cumulative cutoff", {
  uocquantile <- getFromNamespace("uocquantile", "crs")
  x <- ordered(c("low", "low", "mid", "high"), levels = c("low", "mid", "high"))
  expect_identical(as.character(uocquantile(x, 0.5)), "low")
  expect_identical(as.character(uocquantile(x, 0.51)), "mid")
  expect_identical(as.character(uocquantile(x, 0.99)), "high")
})

test_that("uocquantile factor branch returns modal level", {
  uocquantile <- getFromNamespace("uocquantile", "crs")
  x <- factor(c("b", "a", "b", "c"), levels = c("a", "b", "c"))
  expect_identical(as.character(uocquantile(x, 0.5)), "b")
})

test_that("uocquantile numeric branch delegates to quantile", {
  uocquantile <- getFromNamespace("uocquantile", "crs")
  x <- c(1, 2, 4, 8)
  expect_equal(unname(uocquantile(x, 0.25)), unname(stats::quantile(x, probs = 0.25)))
})
