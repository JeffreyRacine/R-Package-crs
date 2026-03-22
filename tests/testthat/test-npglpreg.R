test_that("npglpreg basic estimation works", {
  set.seed(42)
  n <- 100
  x <- runif(n)
  z <- rbinom(n, 1, .5)
  y <- cos(2*pi*x) + z + rnorm(n, sd=.5)
  z <- factor(z)
  
  # Basic formula usage with no CV to be fast
  # We need to provide bws and degree if cv="none"
  # Or we can let it run CV but check that it's silent.
  
  # Let's try a quick run with cv="none" first to check basics
  model <- npglpreg(y ~ x + z, 
                    cv="none", 
                    bws=c(0.1, 0.5), 
                    degree=c(1), 
                    display.warnings=FALSE, 
                    display.nomad.progress=FALSE)
  
  expect_s3_class(model, "npglpreg")
  expect_equal(length(fitted(model)), n)
})

test_that("npglpreg output suppression works", {
  set.seed(42)
  n <- 50
  x <- runif(n)
  y <- x + rnorm(n, sd=0.1)
  
  # Run with CV but suppress output
  # We expect no output to console. 
  # Capturing output is tricky in testthat sometimes but we can check if it runs without error.
  
  out <- capture.output(
    model <- npglpreg(y ~ x, 
                      cv="bandwidth", 
                      degree=c(1),
                      nmulti=1,
                      max.bb.eval=10, # Short run
                      display.warnings=FALSE, 
                      display.nomad.progress=FALSE)
  )
  if(length(out) > 0) {
      print("Captured output:")
      print(out)
  }
  expect_equal(length(out), 0)
  expect_s3_class(model, "npglpreg")
})

test_that("npglpreg prediction and summary works", {
  set.seed(42)
  n <- 50
  x <- runif(n)
  y <- x^2 + rnorm(n, sd=0.1)
  
  model <- npglpreg(y ~ x, 
                    cv="none", 
                    bws=0.1, 
                    degree=2, 
                    display.warnings=FALSE, 
                    display.nomad.progress=FALSE)
  
  # Summary
  expect_output(summary(model))
  
  # Prediction
  newdata <- data.frame(x = seq(0, 1, length=10))
  fit <- predict(model, newdata=newdata)
  expect_equal(length(fit), 10)
})

test_that("npglpreg degree-bandwidth CV works", {
  set.seed(42)
  n <- 30 # Small n for speed
  x <- runif(n)
  y <- sin(2*pi*x) + rnorm(n, sd=0.1)
  
  # Full CV
  model <- npglpreg(y ~ x, 
                    cv="degree-bandwidth", 
                    nmulti=1,
                    max.bb.eval=50, # Sufficient for small n
                    display.warnings=FALSE, 
                    display.nomad.progress=FALSE)
  
  expect_s3_class(model, "npglpreg")
  expect_true(model$degree >= 0)
  expect_true(all(model$bws > 0))
})

test_that("npglpreg NOMAD nmulti stores true outer restart results", {
  set.seed(7)
  n <- 30
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)

  model <- npglpreg(
    y ~ x,
    cv = "degree-bandwidth",
    nmulti = 3,
    max.bb.eval = 20,
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  expect_s3_class(model, "npglpreg")
  expect_identical(model$nmulti, 3)
  expect_identical(model$nomad.nmulti, 3)
  expect_equal(nrow(model$nomad.starts), 3L)
  expect_length(model$nomad.restart.results, 3L)
  expect_length(model$nomad.restart.fval, 3L)
  expect_true(all(is.finite(model$nomad.restart.fval)))
  expect_identical(as.numeric(model$nomad.restart.fval),
                   as.numeric(vapply(model$nomad.restart.results, `[[`, numeric(1L), "objective")))
  expect_identical(model$nomad.best.restart, which.min(model$nomad.restart.fval))
  expect_equal(as.numeric(model$fv), as.numeric(model$nomad.restart.fval[model$nomad.best.restart]))
})

test_that("npglpreg NOMAD progress uses one rich managed line", {
  test_time_values <- function(values) {
    i <- 0L
    force(values)
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  }

  set.seed(42)
  n <- 30L
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)

  traced <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      npglpreg(
        y ~ x,
        cv = "degree-bandwidth",
        nmulti = 3,
        max.bb.eval = 20,
        display.warnings = FALSE,
        display.nomad.progress = TRUE
      )
    ),
    force_renderer = "legacy",
    now = test_time_values(seq(0, 20, by = 0.25)),
    interactive = TRUE
  )

  lines <- vapply(traced$trace, `[[`, character(1L), "line")

  expect_true(any(grepl("^\\[crs\\] Selecting polynomial degree and bw\\.\\.\\.", lines)))
  expect_true(any(grepl("multistart 1/3", lines, fixed = TRUE)))
  expect_true(any(grepl("eval 1", lines, fixed = TRUE)))
  expect_true(any(grepl("deg (", lines, fixed = TRUE)))
  expect_true(any(grepl("best (", lines, fixed = TRUE)))
  expect_true(any(grepl("fv=", lines, fixed = TRUE)))
  expect_false(any(grepl("Calling NOMAD", lines, fixed = TRUE)))
  expect_false(any(grepl("^\\[crs\\] fv = ", lines)))
})

test_that("npglpreg NOMAD progress does not change fixed-seed results", {
  set.seed(19)
  n <- 30L
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  dat <- data.frame(y = y, x = x)

  silent <- suppressWarnings(
    npglpreg(
      y ~ x,
      data = dat,
      cv = "degree-bandwidth",
      nmulti = 2,
      max.bb.eval = 20,
      display.warnings = FALSE,
      display.nomad.progress = FALSE
    )
  )

  noisy <- suppressWarnings(
    npglpreg(
      y ~ x,
      data = dat,
      cv = "degree-bandwidth",
      nmulti = 2,
      max.bb.eval = 20,
      display.warnings = FALSE,
      display.nomad.progress = TRUE
    )
  )

  expect_equal(as.numeric(noisy$degree), as.numeric(silent$degree))
  expect_equal(as.numeric(noisy$bws), as.numeric(silent$bws))
  expect_equal(as.numeric(noisy$fv), as.numeric(silent$fv), tolerance = 1e-12)
  expect_identical(noisy$nomad.best.restart, silent$nomad.best.restart)
  expect_equal(as.numeric(noisy$nomad.restart.fval),
               as.numeric(silent$nomad.restart.fval),
               tolerance = 1e-12)
})

test_that("npglpreg categorical predictors work", {
  set.seed(42)
  n <- 50
  x <- runif(n)
  z <- factor(sample(c("A", "B", "C"), n, replace=TRUE))
  y <- x + as.numeric(z) + rnorm(n, sd=0.1)
  
  model <- npglpreg(y ~ x + z, 
                    cv="none", 
                    bws=c(0.1, 0.05), 
                    degree=1, 
                    display.warnings=FALSE, 
                    display.nomad.progress=FALSE)
  
  expect_s3_class(model, "npglpreg")
  expect_equal(length(model$bws), 2)
  expect_equal(length(fitted(model)), n)
  
  # Check prediction with factors
  newdata <- data.frame(x=0.5, z=factor("B", levels=levels(z)))
  fit <- predict(model, newdata=newdata)
  expect_equal(length(fit), 1)
})

test_that("npglpreg ridging diagnostics honor verbose", {
  set.seed(1)
  n <- 30
  x <- rep(seq(0, 1, length.out = 5), each = 6)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)
  dat <- data.frame(y = y, x = x)

  quiet <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      npglpreg(
        y ~ x,
        data = dat,
        cv = "bandwidth",
        degree = 4,
        nmulti = 1,
        max.bb.eval = 20,
        display.nomad.progress = TRUE,
        display.warnings = TRUE,
        verbose = FALSE
      )
    ),
    force_renderer = "legacy",
    interactive = TRUE
  )

  verbose <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      npglpreg(
        y ~ x,
        data = dat,
        cv = "bandwidth",
        degree = 4,
        nmulti = 1,
        max.bb.eval = 20,
        display.nomad.progress = TRUE,
        display.warnings = TRUE,
        verbose = TRUE
      )
    ),
    force_renderer = "legacy",
    interactive = TRUE
  )

  quiet_lines <- vapply(quiet$trace, `[[`, character(1L), "line")
  verbose_lines <- vapply(verbose$trace, `[[`, character(1L), "line")

  expect_s3_class(quiet$value, "npglpreg")
  expect_s3_class(verbose$value, "npglpreg")
  expect_false(any(grepl("ridging required for inversion", quiet_lines, fixed = TRUE)))
  expect_true(any(grepl("ridging required for inversion", verbose_lines, fixed = TRUE)))
})

test_that("plot.npglpreg bootstrap progress stays visible", {
  test_time_values <- function(values) {
    i <- 0L
    force(values)
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  }
  set.seed(9)
  old_opts <- options(
    crs.progress.start.grace.known.sec = 0,
    crs.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + x2 + rnorm(n, sd = 0.1)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  model <- npglpreg(
    y ~ x1 + x2,
    data = dat,
    cv = "none",
    bws = c(0.2, 0.2),
    degree = c(1, 1),
    display.warnings = FALSE,
    display.nomad.progress = FALSE
  )

  traced <- capture_crs_progress_shadow_trace(
    suppressWarnings(
      plot(
        model,
        mean = TRUE,
        ci = TRUE,
        plot.errors.boot.num = 3,
        plot.behavior = "data",
        num.eval = 8,
        display.warnings = FALSE,
        display.nomad.progress = TRUE
      )
    ),
    force_renderer = "legacy",
    now = test_time_values(seq(0, 20, by = 0.25)),
    interactive = TRUE
  )

  lines <- vapply(traced$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("Plot bootstrap", lines, fixed = TRUE)))
  expect_true(any(grepl("3/3", lines, fixed = TRUE)))
  expect_type(traced$value, "list")
})
