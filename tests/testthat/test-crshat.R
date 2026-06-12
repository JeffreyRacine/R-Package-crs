test_that("crshat matrix and apply outputs reproduce CRS predictions", {
  set.seed(6101)
  d <- data.frame(
    x1 = runif(34),
    x2 = runif(34),
    z = factor(sample(letters[1:3], 34, TRUE))
  )
  d$y <- sin(2 * pi * d$x1) + d$x2 + as.numeric(d$z) / 5 +
    rnorm(34, sd = 0.03)
  nd <- d[c(2, 7, 12, 18, 24), c("x1", "x2", "z")]

  cases <- list(
    additive = crs(
      y ~ x1 + x2 + z, data = d, cv = "none",
      kernel = FALSE, basis = "additive",
      degree = c(2, 1), segments = c(1, 1),
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    tensor_weighted = crs(
      y ~ x1 + x2, data = d, cv = "none",
      kernel = FALSE, basis = "tensor",
      degree = c(2, 2), segments = c(1, 1),
      weights = seq(0.7, 1.3, length.out = nrow(d)),
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    glp_weighted = crs(
      y ~ x1 + x2, data = d, cv = "none",
      kernel = FALSE, basis = "glp",
      degree = c(2, 2), segments = c(1, 1),
      weights = seq(1.3, 0.7, length.out = nrow(d)),
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    kernel_additive = crs(
      y ~ x1 + x2 + z, data = d, cv = "none",
      kernel = TRUE, basis = "additive",
      degree = c(2, 1), segments = c(1, 1),
      lambda = c(0.35),
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    kernel_tensor_weighted = crs(
      y ~ x1 + x2 + z, data = d, cv = "none",
      kernel = TRUE, basis = "tensor",
      degree = c(2, 2), segments = c(1, 1),
      lambda = c(0.45),
      weights = seq(0.8, 1.2, length.out = nrow(d)),
      display.warnings = FALSE, display.nomad.progress = FALSE
    )
  )

  for (nm in names(cases)) {
    model <- cases[[nm]]
    newdata <- nd[, names(model$xz), drop = FALSE]
    pred <- as.vector(predict(model, newdata = newdata))
    H <- crshat(model, newdata = newdata)
    expect_s3_class(H, "crshat")
    expect_equal(as.vector(H %*% model$y), pred,
                 tolerance = 1e-9, info = nm)
    expect_equal(as.vector(crshat(model, newdata = newdata,
                                  output = "apply")),
                 pred, tolerance = 1e-9, info = nm)

    rhs <- cbind(model$y, model$y + seq_along(model$y) / 100)
    expect_equal(crshat(model, newdata = newdata, y = rhs,
                        output = "apply"),
                 H %*% rhs, tolerance = 1e-9, info = nm)
    C <- crshat(model, newdata = newdata, y = model$y,
                output = "constraint")
    expect_equal(matrix(as.numeric(C), nrow = NROW(C)),
                 matrix(as.numeric(t(H) * as.vector(model$y)), nrow = NCOL(H)),
                 tolerance = 1e-9, info = nm)
  }
})

test_that("crshat supports constant fixed-structure fits", {
  set.seed(6102)
  d <- data.frame(x = runif(20), y = rnorm(20))
  model <- crs(
    y ~ x, data = d, cv = "none",
    kernel = FALSE, degree = 0, segments = 1,
    display.warnings = FALSE, display.nomad.progress = FALSE
  )
  nd <- data.frame(x = c(0.1, 0.5, 0.9))
  H <- crshat(model, newdata = nd)
  expect_equal(as.vector(H %*% model$y),
               as.vector(predict(model, newdata = nd)),
               tolerance = 1e-10)
  expect_false(is.infinite(attr(H, "rcond")))
})

test_that("crshat derivative operators reproduce CRS gradient predictions", {
  set.seed(6104)
  n <- 44
  d <- data.frame(
    x1 = runif(n),
    x2 = runif(n),
    z = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  d$y <- sin(2 * pi * d$x1) + d$x2^2 + 0.2 * (d$z == "b") +
    rnorm(n, sd = 0.03)
  nd <- d[c(3, 9, 16, 25, 31), c("x1", "x2", "z")]

  check_derivative_hat <- function(model, newdata, deriv, deriv.index, label) {
    object.deriv <- model
    object.deriv$deriv <- as.integer(deriv)
    pred <- predict(object.deriv, newdata = newdata)
    target <- as.vector(attr(pred, "deriv.mat")[, deriv.index])
    H <- crshat(model, newdata = newdata, deriv = deriv,
                deriv.index = deriv.index)
    expect_equal(as.vector(H %*% model$y), target,
                 tolerance = 1e-9, info = label)
    expect_equal(as.vector(crshat(model, newdata = newdata,
                                  output = "apply", deriv = deriv,
                                  deriv.index = deriv.index)),
                 target, tolerance = 1e-9, info = label)
  }

  cases <- list(
    additive = crs(
      y ~ x1 + x2, data = d, cv = "none", basis = "additive",
      degree = c(3, 3), segments = c(1, 1),
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    tensor_weighted = crs(
      y ~ x1 + x2, data = d, cv = "none", basis = "tensor",
      degree = c(3, 3), segments = c(1, 1),
      weights = seq(0.8, 1.2, length.out = nrow(d)),
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    glp = crs(
      y ~ x1 + x2, data = d, cv = "none", basis = "glp",
      degree = c(3, 3), segments = c(1, 1),
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    kernel_additive = crs(
      y ~ x1 + x2 + z, data = d, cv = "none", kernel = TRUE,
      basis = "additive", degree = c(3, 3), segments = c(1, 1),
      lambda = 0.35,
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    kernel_tensor = crs(
      y ~ x1 + x2 + z, data = d, cv = "none", kernel = TRUE,
      basis = "tensor", degree = c(3, 3), segments = c(1, 1),
      lambda = 0.45,
      display.warnings = FALSE, display.nomad.progress = FALSE
    ),
    kernel_glp = crs(
      y ~ x1 + x2 + z, data = d, cv = "none", kernel = TRUE,
      basis = "glp", degree = c(3, 3), segments = c(1, 1),
      lambda = 0.4,
      display.warnings = FALSE, display.nomad.progress = FALSE
    )
  )

  for (nm in names(cases)) {
    model <- cases[[nm]]
    newdata <- nd[, names(model$xz), drop = FALSE]
    check_derivative_hat(model, newdata, 1L, 1L, paste(nm, "d1x1"))
    check_derivative_hat(model, newdata, 1L, 2L, paste(nm, "d1x2"))
    check_derivative_hat(model, newdata, 2L, 1L, paste(nm, "d2x1"))
  }

  low.degree <- crs(
    y ~ x1, data = d, cv = "none", degree = 1, segments = 1,
    display.warnings = FALSE, display.nomad.progress = FALSE
  )
  H.zero <- crshat(low.degree, newdata = data.frame(x1 = c(0.2, 0.5)),
                  deriv = 2, deriv.index = 1)
  expect_equal(as.vector(H.zero %*% low.degree$y), c(0, 0),
               tolerance = 1e-12)
})

test_that("crshat fails closed for unsupported routes", {
  set.seed(6103)
  d <- data.frame(x = runif(24), y = rnorm(24))
  mean.model <- crs(
    y ~ x, data = d, cv = "none",
    degree = 2, segments = 1,
    display.warnings = FALSE, display.nomad.progress = FALSE
  )
  q.model <- crs(
    y ~ x, data = d, cv = "none",
    degree = 2, segments = 1, tau = 0.5,
    display.warnings = FALSE, display.nomad.progress = FALSE
  )

  expect_error(crshat(q.model), "mean CRS objects")
  expect_error(crshat(mean.model, deriv = 1, deriv.index = 2),
               "deriv.index is invalid")
  expect_error(crshat(mean.model, output = "constraint"),
               "argument 'y' is required")
  expect_error(crshat(mean.model, y = cbind(mean.model$y, mean.model$y),
                      output = "constraint"),
               "one-column")
})
