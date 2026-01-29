## Kernel function used for the categorical variables z
##
## 2025: Optimized to use vector indexing instead of ifelse for performance.

kernel <- function(Z,
                   z,
                   lambda,
                   is.ordered.z=NULL) {

  if(is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda)) stop(" must provide is.ordered.z, Z, z, and lambda")

  if(!is.ordered.z) {
    ## return(ifelse(Z==z,1,lambda))
    res <- rep(lambda, length(Z))
    res[Z == z] <- 1
    return(res)
  } else {
    ## return(ifelse(Z==z,1,lambda^abs(Z-z)))
    res <- lambda^abs(Z-z)
    res[Z == z] <- 1
    return(res)
  }

}

## Product kernel function. Z is a vector/matrix, z is a scalar/vector

prod.kernel <- function(Z,
                        z,
                        lambda,
                        is.ordered.z=NULL,
                        ...,
                        na.rm) {

  Z <- as.matrix(Z)

  if(is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda)) stop(" must provide is.ordered.z, Z, z, and lambda")
  if(length(is.ordered.z) != NCOL(Z)) stop(" is.ordered.z and Z incompatible")

  num.z <- NCOL(Z)

  if(num.z != NROW(z) || num.z != NROW(lambda)) stop(paste(" incompatible dimensions for Z, z, and lambda (",num.z,",",NROW(z),",",NROW(lambda),")",sep=""))

  ## 2025: Inlined kernel logic to avoid function call overhead in loops
  
  ## prodker <- kernel(Z=Z[,1],z=z[1],lambda=lambda[1],is.ordered.z=is.ordered.z[1])
  
  if(!is.ordered.z[1]) {
    prodker <- rep(lambda[1], length(Z[,1]))
    prodker[Z[,1] == z[1]] <- 1
  } else {
    prodker <- lambda[1]^abs(Z[,1]-z[1])
    prodker[Z[,1] == z[1]] <- 1
  }

  if(num.z > 1) {
    for(i in 2:num.z) {
      ## prodker <- prodker * kernel(Z=Z[,i],z=z[i],lambda=lambda[i],is.ordered.z=is.ordered.z[i])
      if(!is.ordered.z[i]) {
        tmp <- rep(lambda[i], length(Z[,i]))
        tmp[Z[,i] == z[i]] <- 1
      } else {
        tmp <- lambda[i]^abs(Z[,i]-z[i])
        tmp[Z[,i] == z[i]] <- 1
      }
      prodker <- prodker * tmp
    }
  }

  return(prodker)
  
}

# function,mean_runtime_original_ns,mean_runtime_fast_ns,mean_speedup,median_speedup,mean_time_reduction_percent
# kernel,36678.8,13783.6,2.661,1.870,62.4
# prod.kernel,119696.6,64880.3,1.845,1.834,45.8

## Faster kernel: same signature & behavior as original
kernel <- function(Z,
                   z,
                   lambda,
                   is.ordered.z = NULL) {
  
  if (is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda))
    stop(" must provide is.ordered.z, Z, z, and lambda")
  
  # Ensure vector behavior is consistent
  Z <- as.vector(Z)
  
  if (!is.ordered.z) {
    res <- rep.int(lambda, length(Z))
    res[Z == z] <- 1
    res
  } else {
    res <- lambda^abs(Z - z)
    res[Z == z] <- 1
    res
  }
}

## Faster product kernel: same signature & behavior as original
prod.kernel <- function(Z,
                        z,
                        lambda,
                        is.ordered.z = NULL,
                        ...,
                        na.rm) {
  
  Z <- as.matrix(Z)
  
  if (is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda))
    stop(" must provide is.ordered.z, Z, z, and lambda")
  if (length(is.ordered.z) != NCOL(Z))
    stop(" is.ordered.z and Z incompatible")
  
  num.z <- NCOL(Z)
  
  ## Follow original API: z and lambda can be column vectors (p x 1)
  ## NROW() is used in the original; keep the same checks
  if (num.z != NROW(z) || num.z != NROW(lambda))
    stop(paste(" incompatible dimensions for Z, z, and lambda (",
               num.z, ",", NROW(z), ",", NROW(lambda), ")", sep = ""))
  
  ## Inline the kernel logic per column (avoids function-call overhead)
  n <- NROW(Z)
  prodker <- rep.int(1, n)
  
  ## Access z/lambda as vectors to mirror original indexing semantics
  z_vec <- as.vector(z)
  lambda_vec <- as.vector(lambda)
  
  for (j in seq_len(num.z)) {
    if (!is.ordered.z[j]) {
      tmp <- rep.int(lambda_vec[j], n)
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    } else {
      tmp <- lambda_vec[j]^abs(Z[, j] - z_vec[j])
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    }
    prodker <- prodker * tmp
  }
  
  prodker
}

