make_kernel_gradient_gram_data <- function(n = 96L, n_eval = 30L,
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

lm_kernel_gradient_reference <- function(y, P.train, P.deriv, basis, weights,
                                         deriv.ind.vec = NULL) {
  P <- P.train
  if(basis=="tensor") {
    model <- lm(y~P-1, weights=weights)
    vcov.model <- vcov(model)
    D <- P.deriv
    beta <- coef(model)
  } else {
    model <- lm(y~P, weights=weights)
    vcov.model <- vcov(model)[-1,-1,drop=FALSE]
    if(basis=="additive") {
      if(is.null(deriv.ind.vec)) {
        stop(" additive derivative columns are missing")
      }
      D <- P.deriv[,deriv.ind.vec,drop=FALSE]
      beta <- coef(model)[-1][deriv.ind.vec]
      vcov.model <- vcov.model[deriv.ind.vec,deriv.ind.vec,drop=FALSE]
    } else {
      D <- P.deriv
      beta <- coef(model)[-1]
    }
  }

  fit <- drop(D %*% beta)
  se <- sqrt(rowSums((D %*% vcov.model) * D))
  list(deriv=fit,
       se=se,
       rank=model$rank,
       df.residual=model$df.residual)
}

test_that("kernel gradient Gram helper matches lm derivative ingredients", {
  data <- make_kernel_gradient_gram_data()
  K <- matrix(c(2, 1, 2, 1), ncol=2L, byrow=TRUE)

  cases <- expand.grid(
    basis = c("additive", "tensor", "glp"),
    weighted = c(FALSE, TRUE),
    mode = c("train", "eval"),
    deriv.index = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    P.train <- crs:::prod.spline(x=data$x, K=K, knots="quantiles",
                                 basis=row$basis,
                                 display.warnings=FALSE)
    if(identical(row$mode, "train")) {
      P.deriv <- crs:::prod.spline(x=data$x, K=K, xeval=data$x,
                                   knots="quantiles", basis=row$basis,
                                   deriv.index=row$deriv.index, deriv=1L,
                                   display.warnings=FALSE)
      z.eval <- data$z
    } else {
      P.deriv <- crs:::prod.spline(x=data$x, K=K, xeval=data$xeval,
                                   knots="quantiles", basis=row$basis,
                                   deriv.index=row$deriv.index, deriv=1L,
                                   display.warnings=FALSE)
      z.eval <- data$zeval
    }

    z.unique <- crs:::uniquecombs(z.eval)
    ind <- attr(z.unique, "index")
    ind.vals <- unique(ind)

    K.additive <- K
    if(identical(row$basis, "additive")) {
      K.additive[,2] <- K[,2]
      K.additive[K[,1] == 0,2] <- 0
      K.additive[,1] <- K[,1]
      K.additive[K[,1] > 0,1] <- K[K[,1] > 0,1] - 1
      deriv.start <- if (row$deriv.index != 1L) {
        sum(K.additive[seq_len(row$deriv.index - 1L), ]) + 1L
      } else {
        1L
      }
      deriv.end <- deriv.start + sum(K.additive[row$deriv.index, ]) - 1L
      deriv.ind.vec <- deriv.start:deriv.end
    } else {
      deriv.ind.vec <- NULL
    }

    for(j in seq_along(ind.vals)) {
      zz <- ind == ind.vals[j]
      L <- crs:::prod.kernel.matrix(
        Z=data$z,
        z=z.unique[ind.vals[j], ],
        lambda=c(0.35, 0.55),
        is.ordered.z=c(FALSE, TRUE)
      )
      if(row$weighted) L <- data$weights * L

      expected <- lm_kernel_gradient_reference(
        y=data$y,
        P.train=P.train,
        P.deriv=P.deriv[zz, , drop=FALSE],
        basis=row$basis,
        weights=L,
        deriv.ind.vec=deriv.ind.vec
      )
      actual <- crs:::.crs_weighted_ls_deriv_lm(
        P.train=P.train,
        y=data$y,
        P.deriv=P.deriv[zz, , drop=FALSE],
        basis=row$basis,
        weights=L,
        deriv.ind.vec=deriv.ind.vec
      )

      expect_false(is.null(actual))
      expect_equal(actual$deriv, expected$deriv,
                   tolerance=1e-10, ignore_attr=TRUE)
      expect_equal(actual$se, expected$se,
                   tolerance=1e-10, ignore_attr=TRUE)
      expect_equal(actual$rank, expected$rank)
      expect_equal(actual$df.residual, expected$df.residual)
    }
  }
})
