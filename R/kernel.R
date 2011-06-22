## Kernel function used for the categorical variables z

kernel <- function(Z,
                   z,
                   lambda,
                   kernel.type=c("nominal","ordinal")) {

  kernel.type=match.arg(kernel.type)

  if(missing(Z) || missing(z) || missing(lambda)) stop(" must provide Z, z, and lambda")

  if(kernel.type=="nominal") {
    return(ifelse(Z==z,1,lambda))
  } else {
    return(ifelse(Z==z,1,lambda^abs(Z-z)))
  }

}

## Product kernel function. Z is a vector/matrix, z is a scalar/vector

prod.kernel <- function(Z,
                        z,
                        lambda,
                        kernel.type=c("nominal","ordinal")) {

  kernel.type=match.arg(kernel.type)

  Z <- as.matrix(Z)

  if(missing(Z) || missing(z) || missing(lambda)) stop(" must provide Z, z, and lambda")

  num.z <- NCOL(Z)

  if(num.z != NROW(z) || num.z != NROW(lambda)) stop(paste(" incompatible dimensions for Z, z, and lambda (",num.z,",",NROW(z),",",NROW(lambda),")",sep=""))
  
  return(apply(sapply(1:num.z, function(j) {kernel(Z=Z[,j],z=z[j],lambda=lambda[j],kernel.type=kernel.type)}),1,prod))

}
