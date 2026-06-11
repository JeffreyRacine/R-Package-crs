test_that("rgl data route preserves payload without loading rgl", {
  rgl.loaded.before <- "rgl" %in% loadedNamespaces()

  set.seed(42)
  n <- 50
  num.eval <- 8
  x1 <- runif(n)
  x2 <- runif(n)
  y <- cos(2 * pi * x1) + sin(2 * pi * x2) + rnorm(n, sd = 0.25)

  model <- crs(
    y ~ x1 + x2,
    degree = c(3, 3),
    segments = c(1, 1),
    cv = "none",
    basis = "tensor",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  payload <- plot(
    model,
    perspective = TRUE, renderer = "rgl",
    output = "data",
    neval = num.eval
  )

  expect_type(payload, "list")
  expect_length(payload, 1L)
  expect_s3_class(payload[[1]], "data.frame")
  expect_equal(nrow(payload[[1]]), num.eval * num.eval)
  expect_named(payload[[1]][1:2], c("x1", "x2"))

  if (!rgl.loaded.before)
    expect_false("rgl" %in% loadedNamespaces())
})

test_that("rgl surface helper restores option and environment state", {
  old.opts <- options(rgl.useNULL = FALSE, rgl.printRglwidget = FALSE)
  on.exit(options(old.opts), add = TRUE)
  old.env <- Sys.getenv("RGL_USE_NULL", unset = NA_character_)
  Sys.setenv(RGL_USE_NULL = "crs-test-sentinel")
  on.exit({
    if (is.na(old.env)) {
      Sys.unsetenv("RGL_USE_NULL")
    } else {
      Sys.setenv(RGL_USE_NULL = old.env)
    }
  }, add = TRUE)

  options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
  Sys.setenv(RGL_USE_NULL = "TRUE")
  skip_if_not_installed("rgl")
  options(rgl.useNULL = FALSE, rgl.printRglwidget = FALSE)
  Sys.setenv(RGL_USE_NULL = "crs-test-sentinel")

  render <- getFromNamespace(".crs_plot_render_surface_rgl", "crs")

  expect_warning(
    capture.output(
      render(
        x = 1:2,
        y = 1:2,
        z = matrix(c(1, 2, 3, 4), 2, 2),
        xlab = "x",
        ylab = "y",
        zlab = "z",
        main = "test",
        col = matrix("white", 2, 2),
        display.warnings = FALSE
      )
    ),
    NA
  )

  expect_identical(getOption("rgl.useNULL"), FALSE)
  expect_identical(getOption("rgl.printRglwidget"), FALSE)
  expect_identical(Sys.getenv("RGL_USE_NULL", unset = NA_character_),
                   "crs-test-sentinel")
})

test_that("rgl plot route completes in rgl null mode", {
  old.opts <- options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
  on.exit(options(old.opts), add = TRUE)
  old.env <- Sys.getenv("RGL_USE_NULL", unset = NA_character_)
  Sys.setenv(RGL_USE_NULL = "TRUE")
  on.exit({
    if (is.na(old.env)) {
      Sys.unsetenv("RGL_USE_NULL")
    } else {
      Sys.setenv(RGL_USE_NULL = old.env)
    }
  }, add = TRUE)

  skip_if_not_installed("rgl")

  set.seed(43)
  n <- 50
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 - x2 + rnorm(n, sd = 0.1)

  model <- crs(
    y ~ x1 + x2,
    degree = c(3, 3),
    segments = c(1, 1),
    cv = "none",
    basis = "tensor",
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  devices.before <- rgl::rgl.dev.list()
  if (is.null(devices.before))
    devices.before <- integer(0L)

  expect_warning(
    capture.output(
      plot(
        model,
        perspective = TRUE, renderer = "rgl",
        neval = 8
      )
    ),
    NA
  )

  devices.after <- rgl::rgl.dev.list()
  if (is.null(devices.after))
    devices.after <- integer(0L)
  expect_equal(sort(unname(devices.after)), sort(unname(devices.before)))
})
