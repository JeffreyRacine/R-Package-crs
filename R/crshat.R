crshat <- function(object, ...) {
  UseMethod("crshat")
}

.crs_hat_constraint_from_matrix <- function(H, y, where) {
  if(is.null(y)) {
    stop(sprintf("argument 'y' is required when output='constraint' in %s",
                 where),
         call. = FALSE)
  }
  y <- as.matrix(y)
  if(NCOL(y) != 1L) {
    stop(sprintf("output='constraint' currently requires vector or one-column 'y' in %s",
                 where),
         call. = FALSE)
  }
  if(NROW(y) != NCOL(H)) {
    stop(sprintf("length of 'y' must match the number of training columns in %s",
                 where),
         call. = FALSE)
  }
  t(H) * as.vector(y)
}

.crs_hat_prepare_newdata <- function(object, newdata) {
  if(is.null(newdata)) {
    return(object$xz)
  }
  Terms <- delete.response(terms(object))
  model.frame(Terms, newdata, xlev = object$xlevels)
}

.crs_hat_split <- function(object, data) {
  if(isTRUE(object$kernel)) {
    splitFrame(data, factor.to.numeric = TRUE)
  } else {
    splitFrame(data)
  }
}

.crs_hat_ginv_from_core <- function(core, p) {
  if(is.null(core$coefficients)) {
    return(NULL)
  }
  rank <- if(identical(core$method, "chol_gram")) p else core$qr$rank
  if(rank < p) {
    return(NULL)
  }
  if(identical(core$method, "chol_gram")) {
    return(chol2inv(core$chol))
  }
  R <- tryCatch(qr.R(core$qr), error = function(e) NULL)
  if(is.null(R)) {
    return(NULL)
  }
  Ginv.pivot <- chol2inv(R)
  pivot <- core$qr$pivot[seq_len(p)]
  if(is.null(pivot)) {
    pivot <- seq_len(p)
  }
  Ginv <- matrix(0, p, p)
  Ginv[pivot, pivot] <- Ginv.pivot
  Ginv
}

.crs_hat_build_design <- function(P.train, P.eval, basis, prune.index = NULL) {
  if(!is.null(prune.index)) {
    P.train <- P.train[, prune.index, drop = FALSE]
    P.eval <- P.eval[, prune.index, drop = FALSE]
  }
  list(X = .crs_weighted_ls_design(P.train, basis),
       X.eval = .crs_weighted_ls_design(P.eval, basis))
}

.crs_hat_additive_deriv_indices <- function(K, deriv.index) {
  K.additive <- K
  K.additive[, 2L] <- K[, 2L]
  K.additive[K[, 1L] == 0, 2L] <- 0
  K.additive[, 1L] <- K[, 1L]
  K.additive[K[, 1L] > 0, 1L] <- K[K[, 1L] > 0, 1L] - 1L
  deriv.start <- if(deriv.index != 1L) {
    sum(K.additive[.crs_index_block(0L, deriv.index - 1L), ]) + 1L
  } else {
    1L
  }
  deriv.end <- deriv.start + sum(K.additive[deriv.index, ]) - 1L
  deriv.start:deriv.end
}

.crs_hat_build_deriv_design <- function(P.train,
                                        P.deriv,
                                        basis,
                                        K,
                                        deriv.index,
                                        prune.index = NULL) {
  X <- .crs_weighted_ls_design(P.train, basis)
  X.eval <- matrix(0, nrow = NROW(P.deriv), ncol = NCOL(X))

  if(identical(basis, "additive")) {
    deriv.ind.vec <- .crs_hat_additive_deriv_indices(K, deriv.index)
    if(!is.null(prune.index)) {
      keep <- prune.index
      P.train <- P.train[, keep, drop = FALSE]
      P.deriv <- P.deriv[, keep, drop = FALSE]
      deriv.ind.vec <- match(deriv.ind.vec, keep)
      deriv.ind.vec <- deriv.ind.vec[!is.na(deriv.ind.vec)]
      X <- .crs_weighted_ls_design(P.train, basis)
      X.eval <- matrix(0, nrow = NROW(P.deriv), ncol = NCOL(X))
    }
    if(length(deriv.ind.vec) > 0L) {
      X.eval[, 1L + deriv.ind.vec] <- P.deriv[, deriv.ind.vec, drop = FALSE]
    }
  } else if(identical(basis, "glp")) {
    if(!is.null(prune.index)) {
      P.train <- P.train[, prune.index, drop = FALSE]
      P.deriv <- P.deriv[, prune.index, drop = FALSE]
      X <- .crs_weighted_ls_design(P.train, basis)
      X.eval <- matrix(0, nrow = NROW(P.deriv), ncol = NCOL(X))
    }
    X.eval[, -1L] <- P.deriv
  } else {
    if(!is.null(prune.index)) {
      P.train <- P.train[, prune.index, drop = FALSE]
      P.deriv <- P.deriv[, prune.index, drop = FALSE]
      X <- .crs_weighted_ls_design(P.train, basis)
    }
    X.eval <- P.deriv
  }

  list(X = X, X.eval = X.eval)
}

.crs_hat_apply_design <- function(X,
                                  X.eval,
                                  Y,
                                  weights = NULL,
                                  rcond.min = 1e-8,
                                  use.svd.fallback = TRUE) {
  weights.work <- if(is.null(weights)) rep(1, NROW(X)) else as.numeric(weights)
  core <- .crs_weighted_ls_core(
    X = X,
    y = rep(0, NROW(X)),
    weights = weights.work,
    rcond.min = rcond.min,
    allow.fallback = TRUE,
    use.svd.fallback = use.svd.fallback
  )
  Ginv <- .crs_hat_ginv_from_core(core, NCOL(X))
  if(is.null(Ginv)) {
    return(NULL)
  }
  Y <- as.matrix(Y)
  if(NROW(Y) != NROW(X)) {
    stop("number of rows in 'y' must match the number of training rows",
         call. = FALSE)
  }
  rhs <- crossprod(X, weights.work * Y)
  out <- X.eval %*% Ginv %*% rhs
  list(value = out,
       method = core$method,
       rcond = core$rcond,
       rank = NCOL(X))
}

.crs_hat_matrix_design <- function(X,
                                   X.eval,
                                   weights = NULL,
                                   rcond.min = 1e-8,
                                   use.svd.fallback = TRUE) {
  weights.work <- if(is.null(weights)) rep(1, NROW(X)) else as.numeric(weights)
  core <- .crs_weighted_ls_core(
    X = X,
    y = rep(0, NROW(X)),
    weights = weights.work,
    rcond.min = rcond.min,
    allow.fallback = TRUE,
    use.svd.fallback = use.svd.fallback
  )
  Ginv <- .crs_hat_ginv_from_core(core, NCOL(X))
  if(is.null(Ginv)) {
    return(NULL)
  }
  H <- X.eval %*% Ginv %*% t(X)
  H <- sweep(H, 2L, weights.work, "*")
  list(H = H,
       method = core$method,
       rcond = core$rcond,
       rank = NCOL(X))
}

.crs_hat_constant_apply <- function(ntrain, neval, Y, weights = NULL) {
  weights.work <- if(is.null(weights)) rep(1, ntrain) else as.numeric(weights)
  denom <- sum(weights.work)
  if(!is.finite(denom) || denom <= 0) {
    return(NULL)
  }
  Y <- as.matrix(Y)
  fit <- matrix(colSums(weights.work * Y) / denom,
                nrow = neval, ncol = NCOL(Y), byrow = TRUE)
  list(value = fit, method = "constant", rcond = NA_real_, rank = 1L)
}

.crs_hat_constant_matrix <- function(ntrain, neval, weights = NULL) {
  weights.work <- if(is.null(weights)) rep(1, ntrain) else as.numeric(weights)
  denom <- sum(weights.work)
  if(!is.finite(denom) || denom <= 0) {
    return(NULL)
  }
  H <- matrix(rep(weights.work / denom, each = neval),
              nrow = neval, ncol = ntrain)
  list(H = H, method = "constant", rcond = NA_real_, rank = 1L)
}

.crs_hat_output_finish <- function(out, vector.output = TRUE) {
  if(is.null(out)) {
    return(NULL)
  }
  if(vector.output && NCOL(out) == 1L) {
    return(as.vector(out))
  }
  out
}

.crs_hat_min_rcond <- function(rconds) {
  rcond <- suppressWarnings(min(rconds, na.rm = TRUE))
  if(!is.finite(rcond)) NA_real_ else rcond
}

.crs_hat_nonkernel <- function(object,
                               eval.data,
                               y,
                               output,
                               deriv,
                               deriv.index,
                               rcond.min,
                               use.svd.fallback) {
  train.split <- .crs_hat_split(object, object$xz)
  eval.split <- .crs_hat_split(object, eval.data)
  x <- as.matrix(train.split$x)
  xeval <- as.matrix(eval.split$x)
  z <- train.split$z
  zeval <- eval.split$z
  K <- object$K
  basis <- object$basis
  weights <- object$weights
  prune.index <- if(isTRUE(object$prune)) object$prune.index else NULL

  deriv <- as.integer(deriv)
  deriv.index <- as.integer(deriv.index)

  if(any(K[, 1L] > 0) || any(object$include > 0)) {
    P.train <- prod.spline(x = x, z = z, K = K, I = object$include,
                           knots = object$knots, basis = basis,
                           display.warnings = FALSE)
    design <- if(deriv > 0L) {
      if(deriv.index < 1L || deriv.index > NROW(K)) {
        stop("deriv.index is invalid", call. = FALSE)
      }
      if(K[deriv.index, 1L] == 0L || deriv > K[deriv.index, 1L]) {
        X <- .crs_weighted_ls_design(
          if(is.null(prune.index)) P.train else P.train[, prune.index, drop = FALSE],
          basis
        )
        list(X = X, X.eval = matrix(0, nrow = NROW(xeval), ncol = NCOL(X)))
      } else {
        P.deriv <- prod.spline(x = x, z = z, K = K, I = object$include,
                               xeval = xeval, zeval = zeval,
                               knots = object$knots, basis = basis,
                               deriv.index = deriv.index, deriv = deriv,
                               display.warnings = FALSE)
        .crs_hat_build_deriv_design(P.train, P.deriv, basis, K,
                                    deriv.index, prune.index)
      }
    } else {
      P.eval <- prod.spline(x = x, z = z, K = K, I = object$include,
                            xeval = xeval, zeval = zeval,
                            knots = object$knots, basis = basis,
                            display.warnings = FALSE)
      .crs_hat_build_design(P.train, P.eval, basis, prune.index)
    }
    if(identical(output, "apply")) {
      return(.crs_hat_apply_design(design$X, design$X.eval, y, weights,
                                   rcond.min, use.svd.fallback))
    }
    return(.crs_hat_matrix_design(design$X, design$X.eval, weights,
                                  rcond.min, use.svd.fallback))
  }

  if(deriv > 0L) {
    if(identical(output, "apply")) {
      return(list(value = matrix(0, nrow = NROW(xeval),
                                 ncol = NCOL(as.matrix(y))),
                  method = "constant_derivative",
                  rcond = NA_real_, rank = 1L))
    }
    return(list(H = matrix(0, nrow = NROW(xeval), ncol = NROW(x)),
                method = "constant_derivative",
                rcond = NA_real_, rank = 1L))
  }

  if(identical(output, "apply")) {
    return(.crs_hat_constant_apply(ntrain = NROW(x), neval = NROW(xeval),
                                   Y = y, weights = weights))
  }
  .crs_hat_constant_matrix(ntrain = NROW(x), neval = NROW(xeval),
                           weights = weights)
}

.crs_hat_kernel <- function(object,
                            eval.data,
                            y,
                            output,
                            deriv,
                            deriv.index,
                            rcond.min,
                            use.svd.fallback) {
  train.split <- .crs_hat_split(object, object$xz)
  eval.split <- .crs_hat_split(object, eval.data)
  x <- as.matrix(train.split$x)
  xeval <- as.matrix(eval.split$x)
  z <- train.split$z
  zeval <- eval.split$z
  if(!is.null(z)) z <- as.matrix(z)
  if(!is.null(zeval)) zeval <- as.matrix(zeval)
  K <- object$K
  basis <- object$basis
  weights <- object$weights
  ntrain <- NROW(x)
  neval <- NROW(xeval)
  deriv <- as.integer(deriv)
  deriv.index <- as.integer(deriv.index)

  if(is.null(z)) {
    if(any(K[, 1L] > 0)) {
      P.train <- prod.spline(x = x, K = K, knots = object$knots,
                             basis = basis, display.warnings = FALSE)
      design <- if(deriv > 0L) {
        if(deriv.index < 1L || deriv.index > NROW(K)) {
          stop("deriv.index is invalid", call. = FALSE)
        }
        if(K[deriv.index, 1L] == 0L || deriv > K[deriv.index, 1L]) {
          X <- .crs_weighted_ls_design(P.train, basis)
          list(X = X, X.eval = matrix(0, nrow = NROW(xeval), ncol = NCOL(X)))
        } else {
          P.deriv <- prod.spline(x = x, K = K, xeval = xeval,
                                 knots = object$knots, basis = basis,
                                 deriv.index = deriv.index, deriv = deriv,
                                 display.warnings = FALSE)
          .crs_hat_build_deriv_design(P.train, P.deriv, basis, K,
                                      deriv.index)
        }
      } else {
        P.eval <- prod.spline(x = x, K = K, xeval = xeval,
                              knots = object$knots, basis = basis,
                              display.warnings = FALSE)
        .crs_hat_build_design(P.train, P.eval, basis)
      }
      if(identical(output, "apply")) {
        return(.crs_hat_apply_design(design$X, design$X.eval, y, weights,
                                     rcond.min, use.svd.fallback))
      }
      return(.crs_hat_matrix_design(design$X, design$X.eval, weights,
                                    rcond.min, use.svd.fallback))
    }
    if(deriv > 0L) {
      if(identical(output, "apply")) {
        return(list(value = matrix(0, nrow = neval,
                                   ncol = NCOL(as.matrix(y))),
                    method = "constant_derivative",
                    rcond = NA_real_, rank = 1L))
      }
      return(list(H = matrix(0, nrow = neval, ncol = ntrain),
                  method = "constant_derivative",
                  rcond = NA_real_, rank = 1L))
    }
    if(identical(output, "apply")) {
      return(.crs_hat_constant_apply(ntrain = ntrain, neval = neval,
                                     Y = y, weights = weights))
    }
    return(.crs_hat_constant_matrix(ntrain = ntrain, neval = neval,
                                    weights = weights))
  }

  if(any(K[, 1L] > 0)) {
    P.train <- prod.spline(x = x, K = K, knots = object$knots,
                           basis = basis, display.warnings = FALSE)
    P.eval <- if(deriv > 0L) {
      if(deriv.index < 1L || deriv.index > NROW(K)) {
        stop("deriv.index is invalid", call. = FALSE)
      }
      if(K[deriv.index, 1L] == 0L || deriv > K[deriv.index, 1L]) {
        NULL
      } else {
        prod.spline(x = x, K = K, xeval = xeval,
                    knots = object$knots, basis = basis,
                    deriv.index = deriv.index, deriv = deriv,
                    display.warnings = FALSE)
      }
    } else {
      prod.spline(x = x, K = K, xeval = xeval,
                  knots = object$knots, basis = basis,
                  display.warnings = FALSE)
    }
  } else {
    P.train <- matrix(1, nrow = ntrain, ncol = 1L)
    P.eval <- if(deriv > 0L) NULL else matrix(1, nrow = neval, ncol = 1L)
    basis <- "tensor"
  }

  zeval.unique <- uniquecombs(as.matrix(zeval))
  ind.zeval <- attr(zeval.unique, "index")
  ind.zeval.vals <- unique(ind.zeval)

  if(identical(output, "apply")) {
    Y <- as.matrix(y)
    out <- matrix(NA_real_, nrow = neval, ncol = NCOL(Y))
    methods <- character(length(ind.zeval.vals))
    rconds <- rep(NA_real_, length(ind.zeval.vals))
    ranks <- integer(length(ind.zeval.vals))
    for(i in seq_along(ind.zeval.vals)) {
      zz <- ind.zeval == ind.zeval.vals[i]
      L <- prod.kernel.matrix(Z = z,
                              z = zeval.unique[ind.zeval.vals[i], ],
                              lambda = object$lambda,
                              is.ordered.z = object$is.ordered.z)
      if(!is.null(weights)) L <- weights * L
      design <- if(deriv > 0L) {
        if(is.null(P.eval)) {
          X <- .crs_weighted_ls_design(P.train, basis)
          list(X = X, X.eval = matrix(0, nrow = sum(zz), ncol = NCOL(X)))
        } else {
          .crs_hat_build_deriv_design(P.train, P.eval[zz, , drop = FALSE],
                                      basis, K, deriv.index)
        }
      } else {
        .crs_hat_build_design(P.train, P.eval[zz, , drop = FALSE], basis)
      }
      tmp <- .crs_hat_apply_design(design$X, design$X.eval, Y, L,
                                   rcond.min, use.svd.fallback)
      if(is.null(tmp)) return(NULL)
      out[zz, ] <- tmp$value
      methods[i] <- tmp$method
      rconds[i] <- tmp$rcond
      ranks[i] <- tmp$rank
    }
    return(list(value = out,
                method = paste(unique(methods), collapse = "+"),
                rcond = .crs_hat_min_rcond(rconds),
                rank = max(ranks)))
  }

  H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
  methods <- character(length(ind.zeval.vals))
  rconds <- rep(NA_real_, length(ind.zeval.vals))
  ranks <- integer(length(ind.zeval.vals))
  for(i in seq_along(ind.zeval.vals)) {
    zz <- ind.zeval == ind.zeval.vals[i]
    L <- prod.kernel.matrix(Z = z,
                            z = zeval.unique[ind.zeval.vals[i], ],
                            lambda = object$lambda,
                            is.ordered.z = object$is.ordered.z)
    if(!is.null(weights)) L <- weights * L
    design <- if(deriv > 0L) {
      if(is.null(P.eval)) {
        X <- .crs_weighted_ls_design(P.train, basis)
        list(X = X, X.eval = matrix(0, nrow = sum(zz), ncol = NCOL(X)))
      } else {
        .crs_hat_build_deriv_design(P.train, P.eval[zz, , drop = FALSE],
                                    basis, K, deriv.index)
      }
    } else {
      .crs_hat_build_design(P.train, P.eval[zz, , drop = FALSE], basis)
    }
    tmp <- .crs_hat_matrix_design(design$X, design$X.eval, L,
                                  rcond.min, use.svd.fallback)
    if(is.null(tmp)) return(NULL)
    H[zz, ] <- tmp$H
    methods[i] <- tmp$method
    rconds[i] <- tmp$rcond
    ranks[i] <- tmp$rank
  }
  list(H = H,
       method = paste(unique(methods), collapse = "+"),
       rcond = .crs_hat_min_rcond(rconds),
       rank = max(ranks))
}

crshat.crs <- function(object,
                       newdata = NULL,
                       y = NULL,
                       output = c("matrix", "apply", "constraint"),
                       deriv = 0,
                       deriv.index = 1,
                       rcond.min = 1e-8,
                       use.svd.fallback = TRUE,
                       ...) {
  output <- match.arg(output)
  if(!inherits(object, "crs")) {
    stop("object must inherit from class 'crs'", call. = FALSE)
  }
  if(!is.null(object$tau)) {
    stop("crshat currently supports mean CRS objects only (tau must be NULL)",
         call. = FALSE)
  }
  deriv <- as.integer(deriv)
  deriv.index <- as.integer(deriv.index)
  if(length(deriv) != 1L || is.na(deriv) || deriv < 0L) {
    stop("deriv must be a non-negative integer", call. = FALSE)
  }
  if(length(deriv.index) != 1L || is.na(deriv.index) || deriv.index < 1L) {
    stop("deriv.index must be a positive integer", call. = FALSE)
  }
  if(is.null(object$xz) || is.null(object$y)) {
    stop("crshat requires a CRS object carrying training data",
         call. = FALSE)
  }
  if(identical(output, "apply") && is.null(y)) {
    y <- object$y
  }
  if(identical(output, "constraint") && is.null(y)) {
    stop("argument 'y' is required when output='constraint' in crshat",
         call. = FALSE)
  }
  eval.data <- .crs_hat_prepare_newdata(object, newdata)
  fit <- if(isTRUE(object$kernel)) {
    .crs_hat_kernel(object, eval.data, y, output, deriv, deriv.index, rcond.min,
                    use.svd.fallback)
  } else {
    .crs_hat_nonkernel(object, eval.data, y, output, deriv, deriv.index, rcond.min,
                       use.svd.fallback)
  }
  if(is.null(fit)) {
    stop("failed to construct CRS hat operator for this fitted object",
         call. = FALSE)
  }

  if(identical(output, "apply")) {
    return(.crs_hat_output_finish(fit$value))
  }

  H <- fit$H
  class(H) <- c("crshat", "matrix")
  attr(H, "object.class") <- class(object)
  attr(H, "train.n") <- NROW(object$xz)
  attr(H, "eval.n") <- NROW(eval.data)
  attr(H, "kernel") <- object$kernel
  attr(H, "basis") <- object$basis
  attr(H, "degree") <- object$degree
  attr(H, "segments") <- object$segments
  attr(H, "lambda") <- object$lambda
  attr(H, "prune") <- object$prune
  attr(H, "prune.index") <- object$prune.index
  attr(H, "deriv") <- deriv
  attr(H, "deriv.index") <- deriv.index
  attr(H, "method") <- fit$method
  attr(H, "rcond") <- fit$rcond
  attr(H, "rank") <- fit$rank
  attr(H, "call") <- match.call(expand.dots = FALSE)
  if(!is.null(y)) {
    Hy <- H %*% as.matrix(y)
    attr(H, "Hy") <- .crs_hat_output_finish(Hy)
  }
  if(identical(output, "constraint")) {
    return(.crs_hat_constraint_from_matrix(H, y, "crshat"))
  }
  H
}

print.crshat <- function(x, ...) {
  cat("\nCRS hat operator:", NROW(x), "evaluation x", NCOL(x), "training\n")
  cat("basis:", attr(x, "basis"),
      "| kernel:", isTRUE(attr(x, "kernel")),
      "| method:", attr(x, "method"), "\n")
  invisible(x)
}
