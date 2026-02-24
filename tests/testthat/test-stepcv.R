test_that("dropterm.default handles empty scope without loop-index artifacts", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)

  res <- crs:::dropterm.default(
    object = fit,
    scope = character(0),
    test = "none",
    display.warnings = FALSE
  )

  expect_s3_class(res, "anova")
  expect_identical(nrow(res), 1L)
  expect_true(is.na(res$Df[1L]))
})

test_that("dropterm.glm handles empty scope without loop-index artifacts", {
  fit <- glm(am ~ wt + hp, data = mtcars, family = binomial())

  res <- crs:::dropterm.glm(
    object = fit,
    scope = character(0),
    test = "none",
    display.warnings = FALSE
  )

  expect_s3_class(res, "anova")
  expect_identical(nrow(res), 1L)
  expect_identical(rownames(res), "<none>")
  expect_true(is.na(res$Df[1L]))
})

test_that("stepCV update path evaluates call objects robustly", {
  fit <- lm(mpg ~ wt, data = mtcars)

  out <- crs:::stepCV(
    object = fit,
    scope = list(lower = mpg ~ 1, upper = mpg ~ wt + hp),
    direction = "both",
    trace = 0,
    steps = 1,
    display.warnings = FALSE
  )

  expect_s3_class(out, "lm")
  expect_true("anova" %in% names(out))
})
