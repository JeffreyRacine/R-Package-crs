make_kernel_predict_gram_data <- function(n = 96L, n_eval = 30L,
                                          seed = 20260610L) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- sample(1:3, n, replace = TRUE)
  z2 <- sample(1:4, n, replace = TRUE)
  y <- sin(2 * pi * x1) + 0.5 * cos(2 * pi * x2) +
    0.2 * z1 + 0.15 * z2 + rnorm(n, sd = 0.10)
  eval_idx <- seq(2L, n, length.out = n_eval)
  list(
    x = cbind(x1 = x1, x2 = x2),
    z = cbind(z1 = z1, z2 = z2),
    y = y,
    weights = seq(0.7, 1.3, length.out = n),
    xeval = cbind(x1 = x1[eval_idx], x2 = x2[eval_idx]),
    zeval = cbind(z1 = z1[eval_idx], z2 = z2[eval_idx])
  )
}

lm_kernel_predict_reference <- function(y, P.train, P.eval, basis, weights,
                                        hat.rows = NULL) {
  P <- P.train
  if(basis=="additive" || basis=="glp") {
    model <- lm(y~P, weights=weights)
  } else {
    model <- lm(y~P-1, weights=weights)
  }
  P <- P.eval
  pred <- predict(model, newdata=data.frame(as.matrix(P)),
                  interval="confidence", se.fit=TRUE)
  out <- cbind(pred[[1]], se=pred[[2]])
  hat <- NULL
  if(!is.null(hat.rows)) hat <- hatvalues(model)[hat.rows]
  list(fitted.values=out,
       hatvalues=hat,
       rank=model$rank,
       df.residual=model$df.residual)
}

test_that("kernel prediction Gram helper matches lm prediction ingredients", {
  data <- make_kernel_predict_gram_data()
  K <- matrix(c(2, 1, 2, 1), ncol=2L, byrow=TRUE)

  cases <- expand.grid(
    basis = c("additive", "tensor", "glp"),
    weighted = c(FALSE, TRUE),
    mode = c("train", "eval"),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    P.train <- crs:::prod.spline(x=data$x, K=K, knots="quantiles",
                                 basis=row$basis,
                                 display.warnings=FALSE)
    if(identical(row$mode, "train")) {
      P.eval <- P.train
      z.eval <- data$z
    } else {
      P.eval <- crs:::prod.spline(x=data$x, K=K, xeval=data$xeval,
                                  knots="quantiles", basis=row$basis,
                                  display.warnings=FALSE)
      z.eval <- data$zeval
    }

    z.unique <- crs:::uniquecombs(z.eval)
    ind <- attr(z.unique, "index")
    ind.vals <- unique(ind)

    for(j in seq_along(ind.vals)) {
      zz <- ind == ind.vals[j]
      L <- crs:::prod.kernel.matrix(
        Z=data$z,
        z=z.unique[ind.vals[j], ],
        lambda=c(0.35, 0.55),
        is.ordered.z=c(FALSE, TRUE)
      )
      if(row$weighted) L <- data$weights * L

      expected <- lm_kernel_predict_reference(
        y=data$y,
        P.train=P.train,
        P.eval=P.eval[zz, , drop=FALSE],
        basis=row$basis,
        weights=L,
        hat.rows=if(identical(row$mode, "train")) which(zz) else NULL
      )
      actual <- crs:::.crs_weighted_ls_predict_lm(
        P.train=P.train,
        y=data$y,
        P.eval=P.eval[zz, , drop=FALSE],
        basis=row$basis,
        weights=L,
        hat.rows=if(identical(row$mode, "train")) which(zz) else NULL
      )

      expect_false(is.null(actual))
      expect_equal(actual$fitted.values, expected$fitted.values,
                   tolerance=1e-10, ignore_attr=TRUE)
      expect_equal(actual$hatvalues, expected$hatvalues,
                   tolerance=1e-10, ignore_attr=TRUE)
      expect_equal(actual$rank, expected$rank)
      expect_equal(actual$df.residual, expected$df.residual)
    }
  }
})
