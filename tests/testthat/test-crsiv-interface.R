iv_interface_fixture <- function(n = 48L) {
  set.seed(20260720)
  v <- rnorm(n, sd = 0.27)
  eps <- rnorm(n, sd = 0.05)
  w <- rnorm(n)
  x <- rnorm(n)
  z <- 0.2 * w + 0.15 * x + v
  y <- z^2 + 0.25 * x - 0.5 * v + eps
  data.frame(y = y, z = z, w = w, x = x,
             wt = seq(0.8, 1.2, length.out = n))
}

iv_interface_args <- function() {
  list(
    cv = "none",
    basis = "additive",
    degree = 2,
    segments = 2,
    nmulti = 1,
    display.nomad.progress = FALSE,
    display.warnings = FALSE
  )
}

test_that("IV formula parser is structural and collision safe", {
  parse_iv <- getFromNamespace(".crs_iv_parse_formula", "crs")
  training_frame <- getFromNamespace(".crs_iv_training_frame", "crs")
  compile_training <- getFromNamespace(".crs_iv_compile_training", "crs")

  parsed <- parse_iv(y ~ I(z^2) | w | x, "probe()")
  expect_identical(parsed$roles$z$labels, "I(z^2)")
  expect_identical(parsed$roles$w$labels, "w")
  expect_identical(parsed$roles$x$labels, "x")

  d <- iv_interface_fixture(12L)
  wrapper <- function(formula, data) {
    mc <- match.call(expand.dots = FALSE)
    parsed <- parse_iv(formula, "probe()")
    training <- training_frame(
      parsed, mc, data, "probe()", parent.frame(), list()
    )
    compile_training(parsed, training, "probe()")
  }
  repeated <- wrapper(y ~ z | z | z, d)
  expect_identical(names(repeated$roles$z), "crsiv_z_001")
  expect_identical(names(repeated$roles$w), "crsiv_w_001")
  expect_identical(names(repeated$roles$x), "crsiv_x_001")
  expect_equal(repeated$roles$z[[1L]], repeated$roles$w[[1L]])

  expect_error(parse_iv(y ~ . | w, "probe()"), "does not support '\\.'")
  expect_error(parse_iv(y ~ z * x | w, "probe()"), "interaction operators")
  expect_error(parse_iv(y ~ offset(z) | w, "probe()"), "offset")
  expect_error(parse_iv(y ~ z, "probe()"), "must be 'y ~ z \\| w'")
  expect_error(parse_iv(y ~ z | w | x | x, "probe()"), "must be 'y ~ z \\| w'")
})

test_that("IV formula compiler supports formula environments and scalar terms", {
  parse_iv <- getFromNamespace(".crs_iv_parse_formula", "crs")
  training_frame <- getFromNamespace(".crs_iv_training_frame", "crs")
  compile_training <- getFromNamespace(".crs_iv_compile_training", "crs")

  wrapper <- function(formula, data = NULL) {
    mc <- match.call(expand.dots = FALSE)
    parsed <- parse_iv(formula, "probe()")
    training <- training_frame(
      parsed, mc, data, "probe()", parent.frame(), list()
    )
    compile_training(parsed, training, "probe()")
  }

  d <- iv_interface_fixture(12L)
  d[["endogenous value"]] <- d$z
  d$group <- factor(rep(c("a", "b"), 6L))
  d$ordered_group <- ordered(rep(c("low", "high"), 6L),
                             levels = c("low", "high"))
  compiled <- wrapper(
    y ~ I(`endogenous value`^2) | group | ordered_group,
    d
  )
  expect_equal(compiled$roles$z[[1L]], d[["endogenous value"]]^2,
               tolerance = 0)
  expect_true(is.factor(compiled$roles$w[[1L]]))
  expect_true(is.ordered(compiled$roles$x[[1L]]))

  environment_fit <- local({
    y <- d$y
    z <- d$z
    w <- d$w
    wrapper(y ~ z | w)
  })
  expect_equal(environment_fit$response, d$y, tolerance = 0)
  expect_equal(environment_fit$roles$z[[1L]], d$z, tolerance = 0)

  expect_error(
    wrapper(y ~ stats::poly(z, 2) | w, d),
    "produces a matrix"
  )
})

test_that("formula data validation cannot be masked by ambient variables", {
  d <- iv_interface_fixture(20L)
  z <- d$z
  malformed <- d[c("y", "w")]
  expect_error(
    crsiv(y ~ z | w, data = malformed),
    "missing required variable.*'z'"
  )
})

test_that("crsiv formula and native routes are numerically identical", {
  d <- iv_interface_fixture()
  args <- iv_interface_args()
  native <- do.call(
    crsiv,
    c(list(y = d$y, z = d["z"], w = d["w"],
           method = "Landweber-Fridman", iterate.max = 3L), args)
  )
  formula <- do.call(
    crsiv,
    c(list(y = y ~ z | w, data = d,
           method = "Landweber-Fridman", iterate.max = 3L), args)
  )

  expect_equal(unname(formula$phi), unname(native$phi), tolerance = 0)
  expect_equal(unname(formula$phi.mat), unname(native$phi.mat), tolerance = 0)
  expect_equal(formula$norm.stop, native$norm.stop, tolerance = 0)
  expect_identical(formula$num.iterations, native$num.iterations)
  expect_identical(formula$convergence, native$convergence)
  expect_s3_class(formula, "crsiv")
  expect_identical(formula$iv$schema_version, 1L)
  expect_identical(formula$iv$roles$z$labels, "z")
  expect_identical(formula$iv$roles$w$labels, "w")
  expect_true(formula$iv$evaluation$is_training_grid)
  expect_identical(formula$iv$public_call[[1L]], quote(crsiv))
  expect_identical(deparse1(formula$iv$public_call[[2L]]), "y ~ z | w")
  expect_identical(environment(formula$iv$public_call[[2L]]), emptyenv())
  expect_identical(native$iv$public_call[[1L]], quote(crsiv))
})

test_that("crsiv formula preserves Tikhonov, quantile, and mixed-factor routes", {
  d <- iv_interface_fixture()
  fixed <- iv_interface_args()

  tikh.native <- do.call(
    crsiv,
    c(list(y = d$y, z = d["z"], w = d["w"],
           method = "Tikhonov", alpha = 1), fixed)
  )
  tikh.formula <- do.call(
    crsiv,
    c(list(y = y ~ z | w, data = d,
           method = "Tikhonov", alpha = 1), fixed)
  )
  expect_equal(unname(tikh.formula$phi), unname(tikh.native$phi),
               tolerance = 0)
  expect_identical(tikh.formula$alpha, tikh.native$alpha)

  quantile.native <- do.call(
    crsiv,
    c(list(y = d$y, z = d["z"], w = d["w"], tau = 0.5,
           method = "Landweber-Fridman", iterate.max = 2L), fixed)
  )
  quantile.formula <- do.call(
    crsiv,
    c(list(y = y ~ z | w, data = d, tau = 0.5,
           method = "Landweber-Fridman", iterate.max = 2L), fixed)
  )
  expect_equal(unname(quantile.formula$phi),
               unname(quantile.native$phi), tolerance = 0)
  expect_equal(quantile.formula$norm.stop,
               quantile.native$norm.stop, tolerance = 0)

  d$instrument_group <- factor(ifelse(seq_len(nrow(d)) %% 2L,
                                      "odd", "even"))
  mixed <- fixed
  mixed$lambda <- 0.3
  factor.native <- do.call(
    crsiv,
    c(list(y = d$y, z = d["z"],
           w = d[c("w", "instrument_group")],
           method = "Landweber-Fridman", iterate.max = 2L), mixed)
  )
  factor.formula <- do.call(
    crsiv,
    c(list(y = y ~ z | w + instrument_group, data = d,
           method = "Landweber-Fridman", iterate.max = 2L), mixed)
  )
  expect_equal(unname(factor.formula$phi), unname(factor.native$phi),
               tolerance = 0)
  expect_equal(factor.formula$norm.stop, factor.native$norm.stop,
               tolerance = 0)
  expect_identical(factor.formula$iv$roles$w$classes[[2L]], "factor")
})

test_that("crsivderiv formula and native selected states are identical", {
  d <- iv_interface_fixture()
  args <- iv_interface_args()
  native <- do.call(
    crsivderiv,
    c(list(y = d$y, z = d["z"], w = d["w"], iterate.max = 3L), args)
  )
  formula <- do.call(
    crsivderiv,
    c(list(y = y ~ z | w, data = d, iterate.max = 3L), args)
  )

  expect_equal(unname(formula$phi), unname(native$phi), tolerance = 0)
  expect_equal(unname(formula$phi.prime), unname(native$phi.prime), tolerance = 0)
  expect_equal(unname(formula$phi.mat), unname(native$phi.mat), tolerance = 0)
  expect_equal(unname(formula$phi.prime.mat), unname(native$phi.prime.mat),
               tolerance = 0)
  expect_equal(formula$norm.stop, native$norm.stop, tolerance = 0)
  expect_identical(formula$num.iterations, native$num.iterations)
  expect_identical(formula$iv$selected_state$iteration,
                   formula$num.iterations)
  expect_identical(formula$iv$selected_state$evaluated_state_count,
                   length(formula$norm.stop))
  expect_identical(formula$iv$selected_state$stopping_value,
                   formula$norm.stop[[formula$num.iterations]])
})

test_that("formula x route preserves the native CRS route", {
  d <- iv_interface_fixture()
  args <- iv_interface_args()
  args$degree <- NULL
  args$segments <- NULL
  native <- do.call(
    crsivderiv,
    c(list(y = d$y, z = d["z"], w = d["w"], x = d["x"],
           iterate.max = 3L), args)
  )
  formula <- do.call(
    crsivderiv,
    c(list(y = y ~ z | w | x, data = d, iterate.max = 3L), args)
  )

  expect_equal(unname(formula$phi), unname(native$phi), tolerance = 0)
  expect_equal(unname(formula$phi.prime), unname(native$phi.prime), tolerance = 0)
  expect_identical(formula$iv$roles$x$labels, "x")
})

test_that("one model frame aligns subset NA weights and starting values", {
  d <- iv_interface_fixture()
  d$y[c(3L, 17L)] <- NA_real_
  d$wt[9L] <- NA_real_
  d$start <- rep(0, nrow(d))
  keep <- seq_len(nrow(d)) %% 5L != 0L
  args <- iv_interface_args()

  fit <- do.call(
    crsiv,
    c(list(y = y ~ z | w, data = d, subset = keep,
           na.action = na.exclude, weights = d$wt,
           starting.values = d$start,
           method = "Landweber-Fridman", iterate.max = 3L), args)
  )

  compact <- d[keep & complete.cases(d[c("y", "z", "w", "wt", "start")]), ]
  native <- do.call(
    crsiv,
    c(list(y = compact$y, z = compact["z"], w = compact["w"],
           weights = compact$wt, starting.values = compact$start,
           method = "Landweber-Fridman", iterate.max = 3L), args)
  )

  expect_equal(unname(fit$phi), unname(native$phi), tolerance = 0)
  expect_equal(length(fitted(fit)), sum(keep))
  expect_equal(length(residuals(fit)), sum(keep))
  expect_equal(sum(is.na(fitted(fit))), 3L)
  expect_equal(sum(is.na(residuals(fit))), 3L)
  expect_equal(fit$iv$rows$n_input, sum(keep))
  expect_equal(length(fit$iv$rows$omitted), 3L)
})

test_that("selected-state accessors and prediction adjacency are coherent", {
  d <- iv_interface_fixture()
  args <- iv_interface_args()
  iv <- do.call(
    crsiv,
    c(list(y = y ~ z | w, data = d, iterate.max = 3L), args)
  )
  ivd <- do.call(
    crsivderiv,
    c(list(y = y ~ z | w, data = d, iterate.max = 3L), args)
  )

  expect_identical(unname(fitted(iv)), unname(iv$phi))
  expect_identical(unname(predict(iv)), unname(iv$phi))
  expect_equal(unname(residuals(iv)), d$y - unname(iv$phi), tolerance = 0)
  expect_identical(unname(fitted(ivd)), unname(ivd$phi))
  expect_identical(unname(predict(ivd)), unname(ivd$phi))
  expect_identical(unname(predict(ivd, deriv = 1L)),
                   unname(ivd$phi.prime))
  expect_equal(unname(residuals(ivd)), d$y - unname(ivd$phi), tolerance = 0)
  expect_error(predict(iv, newdata = d["z"]), "not yet supported")
})

test_that("native post-fit projection remains available", {
  d <- iv_interface_fixture()
  args <- iv_interface_args()
  fit <- do.call(
    crsiv,
    c(list(y = d$y, z = d["z"], w = d["w"], iterate.max = 3L), args)
  )
  newdata <- data.frame(z = seq(min(d$z), max(d$z), length.out = 9L))
  expect_equal(
    unname(predict(fit, newdata = newdata)),
    unname(getFromNamespace("predict.crs", "crs")(fit, newdata = newdata)),
    tolerance = 0
  )
})

test_that("formula evaluation inputs fail before estimator entry", {
  d <- iv_interface_fixture(20L)
  args <- iv_interface_args()
  expect_error(
    do.call(crsiv, c(list(y = y ~ z | w, data = d,
                           zeval = d["z"]), args)),
    "training-row evaluation only"
  )
  expect_error(
    do.call(crsivderiv, c(list(y = y ~ z | w, data = d,
                                newdata = d), args)),
    "training-row evaluation only"
  )
})

test_that("structured IV summaries report the selected state", {
  d <- iv_interface_fixture()
  fit <- do.call(
    crsivderiv,
    c(list(y = y ~ z | w, data = d, iterate.max = 3L), iv_interface_args())
  )
  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary.crsivderiv")
  expect_identical(summary_fit$selected, fit$num.iterations)
  expect_identical(summary_fit$evaluated, length(fit$norm.stop))
  expect_identical(summary_fit$stopping,
                   fit$norm.stop[[fit$num.iterations]])
  output <- capture.output(print(summary_fit))
  expect_true(any(grepl("Selected iteration", output, fixed = TRUE)))
  expect_true(any(grepl("States evaluated", output, fixed = TRUE)))
  expect_true(any(grepl("y ~ z | w", output, fixed = TRUE)))
})

test_that("IV metadata is compact, serializable, and old-object safe", {
  d <- iv_interface_fixture()
  fit <- do.call(
    crsivderiv,
    c(list(y = y ~ z | w, data = d, iterate.max = 3L),
      iv_interface_args())
  )
  stripped <- fit
  stripped[["iv"]] <- NULL
  metadata.bytes <- as.numeric(object.size(fit)) -
    as.numeric(object.size(stripped))
  expect_lt(metadata.bytes, 65536 + 24 * nrow(d))
  expect_identical(environment(fit$iv$public_formula), emptyenv())

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_identical(restored$iv, fit$iv)
  expect_true(isTRUE(all.equal(restored, fit, check.environment = FALSE)))
  expect_identical(unname(fitted(restored)), unname(fit$phi))

  legacy <- stripped
  expect_identical(fitted(legacy), legacy$phi)
  expect_identical(residuals(legacy), legacy$residuals)
  expect_s3_class(summary(legacy), "summary.crsivderiv")
  expect_silent(capture.output(print(legacy)))
})

test_that("IV plot payloads use selected evaluation state and original y", {
  d <- iv_interface_fixture()
  fit <- do.call(
    crsivderiv,
    c(list(y = y ~ z | w, data = d, iterate.max = 3L), iv_interface_args())
  )
  payload <- getFromNamespace(".crs_plot_payload_iv_deriv", "crs")(
    fit, phi = FALSE
  )
  expect_equal(payload$data[[1L]], sort(d$z), tolerance = 0)
  expect_equal(payload$data$fit, fit$phi.prime[order(d$z)], tolerance = 0)

  overlay <- getFromNamespace(".crs_iv_plot_training_data", "crs")(fit)
  expect_equal(overlay$z, d$z, tolerance = 0)
  expect_equal(overlay$y, d$y, tolerance = 0)
})

test_that("off-training native objects refuse residuals and overlays", {
  d <- iv_interface_fixture()
  args <- iv_interface_args()
  fit <- do.call(
    crsivderiv,
    c(list(y = d$y, z = d["z"], w = d["w"],
           zeval = d["z"], weval = d["w"], iterate.max = 3L), args)
  )
  expect_false(fit$iv$evaluation$is_training_grid)
  expect_error(residuals(fit), "training rows")
  expect_error(
    getFromNamespace(".crs_iv_plot_training_data", "crs")(fit),
    "training rows"
  )
})
