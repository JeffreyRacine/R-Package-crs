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
  
  ## 2025: index.swap is equivalent to order() for permutations
  index.swap <- function(index.input)
  {
    return(order(index.input))
  }
  
  ## 2025: Optimized two.dimen with pre-allocation
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
    
    # Pre-filter d1sets if needed
    nc <- NCOL(d1sets)
    use_d1sets <- d1sets
    if (nc > 1 && d2p > 0) {
      keep <- d1sets[, nc] != d2p
      use_d1sets <- d1sets[keep, , drop = FALSE]
    }
    rs <- rowSums(use_d1sets)
    
    # Count rows for pre-allocation
    # Using tabulate for speed (indices are positive integers)
    # rs contains row sums. 
    max_rs <- if(length(rs) > 0) max(rs) else 0
    # Ensure counts vector is large enough to handle max index d1
    counts <- integer(max(max_rs, d1) + 1)
    if (length(rs) > 0) {
      tab <- table(rs)
      counts[as.integer(names(tab))] <- as.integer(tab)
    }
    
    total_rows <- d2
    
    if(d1 > d2) {
      idx1 <- 1:(d1-d2)
      total_rows <- total_rows + d2 * sum(counts[idx1])
    }
    
    idx2 <- 1:d2
    target <- d1 - idx2 + 1
    total_rows <- total_rows + sum(idx2 * counts[target])
    
    # Allocate result
    d2sets <- matrix(0, nrow=total_rows, ncol=nc+1)
    
    # Initialize first d2 rows
    # 1:nc cols are 0, last col is 1:d2
    d2sets[1:d2, nc+1] <- 1:d2
    
    curr_row <- d2 + 1
    
    # Loop 1
    if(d1-d2 > 0){
      for(i in 1:(d1-d2)){
        # d12 update logic from original (tracking theoretical size vs actual?)
        # d12 <- d12+d2*nd1[i] 
        # Note: nd1 is used for nd2 calc later, but here we use actual counts from d1sets
        
        # Extract rows with sum i
        idx <- which(rs == i)
        n_sub <- length(idx)
        
        if (n_sub > 0) {
          end_row <- curr_row + n_sub * d2 - 1
          
          # Replicate rows d2 times
          d2sets[curr_row:end_row, 1:nc] <- use_d1sets[rep(idx, d2), , drop=FALSE]
          # Last column: rep(0:(d2-1), each=n_sub)
          d2sets[curr_row:end_row, nc+1] <- rep(0:(d2-1), each=n_sub)
          
          curr_row <- end_row + 1
        }
      }
    }
    
    # Loop 2
    for(i in 1:d2){
      # d12 <- d12 + (i*nd1[d1-i+1])
      
      target <- d1 - i + 1
      idx <- which(rs == target)
      n_sub <- length(idx)
      
      if (n_sub > 0) {
        end_row <- curr_row + n_sub * i - 1
        
        d2sets[curr_row:end_row, 1:nc] <- use_d1sets[rep(idx, i), , drop=FALSE]
        d2sets[curr_row:end_row, nc+1] <- rep(0:(i-1), each=n_sub)
        
        curr_row <- end_row + 1
      }
    }
    
    # nd2 calculation (unchanged)
    nd2 <- nd1
    if(d1>1){
      for(j in 1:(d1-1)) {
        nd2[j] <- 0
        limit <- max(0, j-d2+1)
        low <- max(1, limit)
        if (j >= low) nd2[j] <- nd2[j] + sum(nd1[low:j])
        if (limit == 0) nd2[j] <- nd2[j] + 1
      }
    }
    if(d2>1) {
      nd2[d1] <- nd1[d1]
      low <- d1-d2+1
      high <- d1-1
      if(high >= low) nd2[d1] <- nd2[d1] + sum(nd1[low:high])
    }
    else {
      nd2[d1] <- nd1[d1]
    }
    
    ret <- list()
    ret$d12 <- total_rows
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
  ## 2025: Faster list creation for order
  z <- z[do.call('order', lapply(ncol(z):1, function(i) z[,i])), ]
  
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

