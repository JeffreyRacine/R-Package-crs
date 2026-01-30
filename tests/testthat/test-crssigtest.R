test_that("crssigtest works", {
  options(crs.messages=FALSE)
  set.seed(42)
  n <- 100
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd=0.1)
  
  model <- crs(y ~ x1 + x2, degree=c(3,3), segments=c(1,1), cv="none", basis="additive", display.warnings=FALSE, display.nomad.progress=FALSE)
  
  # Test with very few bootstrap reps for speed
  sig <- crssigtest(model, boot.num=9, boot=TRUE)
  
  expect_s3_class(sig, "sigtest.crs")
  expect_length(sig$P, 2)
  expect_true(all(sig$P >= 0 & sig$P <= 1))
})

test_that("crssigtest works without bootstrap", {
  options(crs.messages=FALSE)
  set.seed(42)
  n <- 100
  x1 <- runif(n)
  y <- x1 + rnorm(n, sd=0.1)
  
  model <- crs(y ~ x1, degree=3, segments=1, cv="none", basis="additive", display.warnings=FALSE, display.nomad.progress=FALSE)
  
  sig <- crssigtest(model, boot=FALSE)
  
  expect_s3_class(sig, "sigtest.crs")
  expect_length(sig$P.asy, 1)
})
