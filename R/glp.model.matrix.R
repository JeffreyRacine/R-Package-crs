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
  
  k <-length(X)
  dimen <- numeric()
  for(i in 1:k) {
    dimen[i] <- ncol(X[[i]])
  }
  
  index.swap <- function(index.input)
  {
    index.ret <- 1:length(index.input)
    index <- 1:length(index.input)
    for (i in 1:length(index.input))
    {
      index.ret[i] <- index[index.input==i]
    }
    
    return (index.ret)
  }
  
  two.dimen<- function(d1,d2,d2p,nd1,pd12,d1sets)
  {
    if(d2 == 1) {
      ret <- list()
      ret$d12 <- pd12
      ret$nd1 <- nd1
      ret$d2p <- d2
      ret$sets <- cbind(d1sets,0)
      return(ret)
    }
    d12 <- d2
    d2sets <- cbind(matrix(0,nrow=min(NROW(d1sets),d2),ncol=NCOL(d1sets)),seq(1:d2))
    if(d1-d2>0){
      for(i in 1:(d1-d2)){
        d12 <- d12+d2*nd1[i]
        oneSet <- d1sets
        if(NCOL(d1sets) > 1 & d2p > 0)
          oneSet <- d1sets[d1sets[,NCOL(d1sets)]!=d2p,]  #previous d2
        d2sets <- rbind(d2sets,cbind(matrix(oneSet[rowSums(oneSet)==i,],nrow=NROW(oneSet[rowSums(oneSet)==i,]),ncol=NCOL(d1sets))[rep(1:NROW(oneSet[rowSums(oneSet)==i,]),d2),],rep(0:(d2-1),each=nd1[i])))
      }
    }
    
    for(i in 1:d2){
      d12 <- d12 + (i*nd1[d1-i+1])
      if(nd1[d1-i+1]>0){
        oneSet <- d1sets
        if(NCOL(d1sets) > 1 & d2p > 0)
          oneSet <- d1sets[d1sets[,NCOL(d1sets)]!=d2p,]  #previous d2
        d2sets <- rbind(d2sets,cbind(matrix(oneSet[rowSums(oneSet)==d1-i+1,],nrow=NROW(oneSet[rowSums(oneSet)==d1-i+1,]),ncol=NCOL(d1sets))[rep(1:NROW(oneSet[rowSums(oneSet)==d1-i+1,]),i),],rep(0:(i-1), each=nd1[d1-i+1])))
      }
    }
    
    nd2 <- nd1  ## Calculate nd2
    if(d1>1){
      for(j in 1:(d1-1)) {
        nd2[j] <- 0
        for(i in j:max(0,j-d2+1)) {
          if(i > 0) {
            nd2[j] <- nd2[j] + nd1[i]                  
          }
          else {
            nd2[j] <- nd2[j] + 1  ## nd1[0] always 1
          }
        }
      }
    }
    if(d2>1) {
      nd2[d1] <- nd1[d1]
      for(i in (d1-d2+1):(d1-1)) nd2[d1] <- nd2[d1]+nd1[i]
    }
    else {
      nd2[d1] <- nd1[d1]
    }
    ret <- list()
    ret$d12 <- d12
    ret$nd1 <- nd2
    ret$d2p <- d2
    ret$sets <- d2sets
    
    return(ret)
  }
  
  construct.tensor.prod <- function(dimen.input)
  {
    dimen <- sort(dimen.input, decreasing = TRUE)
    
    nd1 <- rep(1,dimen[1])   ## At the beginning,  we have one for [1, 2, 3, ..., dimen[1]]
    nd1[dimen[1]] <- 0       ## nd1 represents the frequency for every element of [1, 2, 3, ..., dimen[1]]
    ncol.bs <- dimen[1]
    sets <- 1:dimen[1]
    dim(sets)<- c(dimen[1],1)
    d2p <- 0
    
    if(k>1) {
      for(i in 2:k) {
        dim.rt <- two.dimen(dimen[1],dimen[i],d2p,nd1,ncol.bs,sets)
        nd1 <- dim.rt$nd1
        ncol.bs <- dim.rt$d12
        sets <- dim.rt$sets
        d2p <- dim.rt$d2p
      }
      ncol.bs <- dim.rt$d12+k-1
      for(i in 2:k)
      {
        oneRow <- rep(0, NCOL(sets))
        oneRow[i-1] <- dimen[i-1]
        sets <- rbind(sets, oneRow)
      }
    }
    
    return(sets[, index.swap(order(dimen.input, decreasing = TRUE)), drop = FALSE])
  }
  
  z <- construct.tensor.prod(dimen)
  ## if we don't want to get the same order as using the expand.grid method, we can comment the following two lines.
  rownames(z) <- 1:NROW(z)
  z <- z[do.call('order', as.list(data.frame(z[, NCOL(z):1]))), ]
  
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

