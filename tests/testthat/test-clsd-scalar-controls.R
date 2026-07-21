test_that("clsd rejects ambiguous scalar route controls explicitly", {
  invalid <- list(c(TRUE, FALSE), logical(), NA, "neither")

  for (value in invalid) {
    expect_error(
      clsd(elastic.max = value),
      "^elastic\\.max must be TRUE or FALSE$"
    )
    expect_error(
      clsd(NOMAD = value),
      "^NOMAD must be TRUE or FALSE$"
    )
  }
})

test_that("clsd preserves valid scalar condition forms and error precedence", {
  valid <- list(TRUE, FALSE, 0, 1)

  for (value in valid) {
    expect_error(
      clsd(elastic.max = value, NOMAD = FALSE),
      "^ You must provide data$"
    )
    expect_error(
      clsd(elastic.max = FALSE, NOMAD = value),
      "^ You must provide data$"
    )
  }
})

test_that("clsd scalar normalization preserves search routing and elastic caps", {
  capture_route <- function(elastic.max, NOMAD) {
    seen <- NULL

    fit <- testthat::with_mocked_bindings(
      clsd(
        x = seq(-1, 1, length.out = 8L),
        degree.max = 9,
        segments.max = 8,
        elastic.max = elastic.max,
        NOMAD = NOMAD,
        n.integrate = 20,
        display.warnings = FALSE,
        display.nomad.progress = FALSE
      ),
      ls.ml = function(...) {
        seen <<- list(...)
        list(beta = rep(-1, 3), degree = 2, segments = 1, fv = 0)
      },
      .package = "crs"
    )
    expect_s3_class(fit, "clsd")

    seen[c("degree.max", "segments.max", "elastic.max", "NOMAD")]
  }

  exhaustive.elastic <- capture_route(TRUE, FALSE)
  expect_identical(exhaustive.elastic$degree.max, 3)
  expect_identical(exhaustive.elastic$segments.max, 3)
  expect_identical(exhaustive.elastic$elastic.max, TRUE)
  expect_identical(exhaustive.elastic$NOMAD, FALSE)

  exhaustive.fixed <- capture_route(FALSE, FALSE)
  expect_identical(exhaustive.fixed$degree.max, 9)
  expect_identical(exhaustive.fixed$segments.max, 8)
  expect_identical(exhaustive.fixed$elastic.max, FALSE)
  expect_identical(exhaustive.fixed$NOMAD, FALSE)

  nomad <- capture_route(TRUE, TRUE)
  expect_identical(nomad$degree.max, 9)
  expect_identical(nomad$segments.max, 8)
  expect_identical(nomad$elastic.max, TRUE)
  expect_identical(nomad$NOMAD, TRUE)
})
