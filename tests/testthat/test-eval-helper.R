test_that(".crs_eval_call matches eval() for update() call objects", {
  data <- transform(mtcars, y = mpg)
  fit <- lm(y ~ cyl, data = data)
  fit_call <- update(fit, . ~ . + wt, evaluate = FALSE)
  env <- environment(formula(fit))

  fit_new <- crs:::.crs_eval_call(fit_call, env)
  fit_ref <- eval(fit_call, envir = env)

  expect_equal(coef(fit_new), coef(fit_ref), tolerance = 1e-12)
  expect_equal(fitted(fit_new), fitted(fit_ref), tolerance = 1e-12)
})

test_that(".crs_eval_call resolves namespace-qualified call heads", {
  env <- new.env(parent = baseenv())
  env$x <- c(1, 2, 3, 4, 5)

  expect_equal(crs:::.crs_eval_call(quote(stats::median(x)), env), 3)
})

test_that("succeedWithResponse returns expected response-availability flags", {
  tt <- terms(y ~ x1 + x2)

  df_has_response <- data.frame(y = 1:4, x1 = 5:8, x2 = 9:12)
  df_missing_response <- data.frame(x1 = 5:8, x2 = 9:12)

  expect_true(crs:::succeedWithResponse(tt, df_has_response))
  expect_false(crs:::succeedWithResponse(tt, df_missing_response))
})
