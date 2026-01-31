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
