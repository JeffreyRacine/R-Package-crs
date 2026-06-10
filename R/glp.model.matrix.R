## (C) Jeffrey S. Racine July 22 2011

## glp.model.matrix is a modified version of the polym() function
## (stats) combined with the tensor.prod.model.matrix function in
## mgcv. The function accepts a vector of degrees and provides a
## generalized polynomial with varying polynomial order. This can be
## more parsimonious than the tensor product model matrix commonly
## found in the spline literature while retaining solid approximation
## capabilities and will be better conditioned than the tensor product
## (see commented illustration below).

## X is a list of model matrices, from which a generalized local
## polynomial model matrix is to be produced (first column ones).

## Note - if this is used for modeling and derivatives are required,
## one would pass in the derivative model matrix for the variable(s)
## whose derivative is required and pass matrices of zeros of
## identical dimension to the originals for all other variables.

glp.model.matrix <- function(X) {

  k <- length(X)
  dimen <- numeric()
  for(i in seq_len(k)) {
    dimen[i] <- ncol(X[[i]])
  }

  dimen.list <- vector("list", k)
  for(i in seq_len(k)) {
    dimen.list[[i]] <- 0:dimen[i]
  }

  z <- do.call("expand.grid", dimen.list)
  s <- rowSums(z)
  z <- z[(s > 0) & (s <= max(dimen)), , drop = FALSE]

  if(!all(dimen == max(dimen))) {
    for(j in seq_along(dimen)) {
      d <- dimen[j]
      if((d < max(dimen)) && (d > 0)) {
        s <- rowSums(z)
        drop <- (s > d) & (z[, j, drop = FALSE] ==
                             matrix(d, nrow(z), 1, byrow = TRUE))
        z <- z[!drop, , drop = FALSE]
      }
    }
  }

  # X is a list of model matrices, from which a tensor product model matrix is to be produced.
  # e.g. ith row is basically X[[1]][i,]%x%X[[2]][i,]%x%X[[3]][i,], but this routine works
  # column-wise, for efficiency, and does work in compiled code.
  m <- length(X)              ## number to row tensor product
  d <- unlist(lapply(X,ncol)) ## dimensions of each X
  n <- nrow(X[[1]])           ## columns in each X
  zn <- NROW(z)
  X <- as.numeric(unlist(X))  ## append X[[i]]s columnwise
  z <- as.integer(unlist(z))  ## append z's columnwise
  res <- numeric(n*zn)     ## storage for result
  .Call(glp_model_tmm,X,z, res,d,m,n, zn)   ## produce product

  return(matrix(res,nrow=n))

}

##> set.seed(42)
##> n <- 1000
##>
##> degree <- 10
##> nbreak <- 2
##>
##> X1 <- gsl.bs(runif(n),degree=degree,nbreak=nbreak)
##> X2 <- gsl.bs(runif(n),degree=degree,nbreak=nbreak)
##> X3 <- gsl.bs(runif(n),degree=degree,nbreak=nbreak)
##>
##> X <- list()
##> X[[1]] <- X1
##> X[[2]] <- X2
##> X[[3]] <- X3
##>
##>   k<-length(X)
##>   dimen.list <- list()
##>   for(i in 1:k) dimen.list[[i]] <- 0:ncol(X[[i]])
##>
##> B.glp <- glp.model.matrix(X)
##> B.tp <- tensor.prod.model.matrix(X)
##>
##> dim(B.glp)
##[1] 1000  285
##> dim(B.tp)
##[1] 1000 1000
##>
##> all.equal(matrix(B.glp),matrix(B.tp))
##[1] "Attributes: < Component 1: Mean relative difference: 2.508772 >"
##[2] "Numeric: lengths (285000, 1000000) differ"
##>
##> rcond(t(B.glp)%*%B.glp)
##[1] 8.338926e-13
##> rcond(t(B.tp)%*%B.tp)
##[1] 5.92409e-21
