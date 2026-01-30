## No Zero Denominator, used in C code for kernel estimation...

## Original, better, best

# NZD <- function(a) {
#   ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
# }

# NZD <- function(a) {
#   if(length(a) == 1) {
#     if(is.na(a)) return(a)
#     if(a < 0) return(min(-.Machine$double.eps, a))
#     return(max(.Machine$double.eps, a))
#   }
#   ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
# }

# Benches 5.5-7.9 times faster (mean/median) produces identical results

NZD <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1) {
    if (a >= 0) {
      if (a < eps) return(eps)
    } else {
      if (a > -eps) return(-eps)
    }
    return(a)
  }
  idx <- which(abs(a) < eps)
  if (length(idx) > 0)
    a[idx] <- ifelse(a[idx] >= 0, eps, -eps)
  a
}

## New function when only positive values are guaranteed, better, best

# NZD_pos <- function(a) {
#   if(length(a) == 1) {
#     if(is.na(a)) return(a)
#     return(max(.Machine$double.eps, a))
#   }
#   pmax(.Machine$double.eps,a)
# }

# Benches 1.3-1.8 times faster (mean/median) produces identical results

NZD_pos <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1)
    return(if (a < eps) eps else a)
  idx <- which(a < eps)
  if (length(idx) > 0)
    a[idx] <- eps
  a
}

## Integration timing:
# Comparison,Mean Speedup,Median Speedup
# Cumulative Integration (fast_cum vs original_cum),1.54x,1.51x
# Summation Integration (fast_sum vs original_sum),1.55x,1.25x

# integrate.trapezoidal <- function(x,y) {
# 
#   ## This function will compute the cumulative integral at each sample
#   ## realization using the Newton-Cotes trapezoidal rule and the
#   ## cumsum function as we need to compute this in a computationally
#   ## efficient manner. It can be used to return the distribution
#   ## function from the density function etc. We test for unsorted data.
# 
#   n <- length(x)
#   if(x.unsorted <- is.unsorted(x)) {
#     rank.x <- rank(x)
#     order.x <- order(x)
#     y <- y[order.x]
#     x <- x[order.x]
#   }
# 
#   int.vec <- numeric(length(x))
#   int.vec[2:n] <- cumsum((x[2:n] - x[2:n-1]) * (y[2:n] + y[2:n-1]) / 2)
# 
#   if(x.unsorted) {
#     return(int.vec[rank.x])
#   } else {
#     return(int.vec)
#   }
#     
# }

integrate.trapezoidal <- function(x, y) {
  if (is.unsorted(x)) {
    idx <- order(x)
    # We use rank via the order index to avoid calling rank() separately
    rnk <- match(seq_along(x), idx)
    x <- x[idx]
    y <- y[idx]
  } else {
    rnk <- NULL
  }
  
  # Optimization: Compute diffs and sums once
  # x[-1] is x[2:n], x[-length(x)] is x[1:n-1]
  n <- length(x)
  res <- c(0, cumsum((x[-1] - x[-n]) * (y[-1] + y[-n]) / 2))
  
  if (!is.null(rnk)) return(res[rnk])
  return(res)
}

# integrate.trapezoidal.sum <- function(x,y) {
# 
#   ## This function will compute the cumulative integral at each sample
#   ## realization using the Newton-Cotes trapezoidal rule and the
#   ## cumsum function as we need to compute this in a computationally
#   ## efficient manner. It can be used to return the distribution
#   ## function from the density function etc. We test for unsorted data.
# 
#   n <- length(x)
#   if(x.unsorted <- is.unsorted(x)) {
#     rank.x <- rank(x)
#     order.x <- order(x)
#     y <- y[order.x]
#     x <- x[order.x]
#   }
# 
#   return(sum((x[2:n] - x[2:n-1]) * (y[2:n] + y[2:n-1]) / 2))
#     
# }

integrate.trapezoidal.sum <- function(x, y) {
  if (is.unsorted(x)) {
    idx <- order(x)
    x <- x[idx]
    y <- y[idx]
  }
  n <- length(x)
  # Efficient vector multiplication
  return(sum((x[-1] - x[-n]) * (y[-1] + y[-n])) / 2)
}

## This function tests for monotone increasing vectors

is.monotone.increasing <- function(x) {
  ## Sorted and last value > first value
  !is.unsorted(x) && x[length(x)] > x[1]
}

## This function tests for the maximum well-conditioned spline degree.

## Note that increasing the number of breaks, other things equal,
## results in a better-conditioned matrix. Hence we ignore nbreak and
## set it to its minimum (2)

check.max.spline.degree <- function(xdat=NULL,degree=NULL,display.warnings=TRUE) {
  
  if(is.null(xdat)) stop(" xdat must be provided")
  if(is.null(degree)) stop(" degree vector must be provided")
  
  xdat <- as.data.frame(xdat)
  
  if(missing(degree)) stop(" degree vector must be provided")
  
  ill.conditioned <- FALSE
  
  xdat.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})
  numeric.index <- which(xdat.numeric==TRUE)  
  num.numeric <- sum(sapply(1:NCOL(xdat),function(i){is.numeric(xdat[,i])})==TRUE)
  d <- numeric(num.numeric)
  
  if(num.numeric > 0) {
    
    for(i in 1:num.numeric) {
      if(degree[i]>0) {
        X <- gsl.bs(xdat[,numeric.index[i]],degree=degree[i],nbreak=2)
        d[i] <- degree[i]
        if(!is.fullrank(X)) {
          for(j in 1:degree[i]) {
            d[i] <- j
            X <- gsl.bs(xdat[,numeric.index[i]],degree=d[i],nbreak=2)
            if(!is.fullrank(X)) {
              d[i] <- j-1
              break()
            }
          }
        }
        if(d[i] < degree[i]) {
          if(display.warnings) warning(paste("\r Predictor ",i," B-spline basis is ill-conditioned beyond degree ",d[i],": see note in ?npglpreg",sep=""),immediate.=TRUE)
          ill.conditioned <- TRUE
        }
      }
    }
    
  }
  
  attr(ill.conditioned, "degree.max.vec") <- d
  return(ill.conditioned)
  
}

succeedWithResponse <- function(tt, frame){
  !any(class(try(eval(expr = attr(tt, "variables"),
                      envir = frame, enclos = NULL), silent = TRUE)) == "try-error")
}

## Utility function to divide explanatory variables into
## factors/numeric, strip off names etc.

# splitFrame <- function(xz, factor.to.numeric=FALSE) {
#   
#   if(missing(xz)) stop(" you must provide xz data")
#   if(!is.data.frame(xz)) stop(" xz must be a data frame")
# 
#   xznames <- names(xz)
#   
#   IND <- logical()
# 
#   for(i in 1:NCOL(xz)) IND[i] <- is.factor(xz[,i])
# 
#   x <- xz[,!IND,drop=FALSE]
#   num.x <- ncol(x)
# 
#   ## We require at least one continuous predictor to conduct spline
#   ## smoothing, but there may/may not be factors.
# 
#   if(num.x == 0) stop(" can't fit spline surfaces with no continuous predictors")
# 
#   xnames <- xznames[!IND]
# 
#   is.ordered.z <- NULL
#   
#   if(any(IND)) {
#     is.ordered.z <- logical()
#     for(i in 1:NCOL(xz[,IND,drop=FALSE])) is.ordered.z[i] <- is.ordered((xz[,IND,drop=FALSE])[,i])
#     if(!factor.to.numeric) {
#       z <- data.frame(xz[,IND,drop=FALSE])
#     } else {
#       ## If factor.to.numeric crudely convert factors to numeric.
#       z <- matrix(NA,NROW(xz),NCOL(xz[,IND,drop=FALSE]))
#       ## To "revert" a factor f to its original numeric values,
#       ## as.numeric(levels(f))[f] is recommended. Problem is that for
#       ## character strings it produces a warning message. No idea how
#       ## to test for this so dropping for the moment. Will affect
#       ## ordered types.
#       for(i in 1:NCOL(xz[,IND,drop=FALSE])) {
#         suppressWarnings(z[,i] <- as.numeric(levels((xz[,IND,drop=FALSE])[,i]))[(xz[,IND,drop=FALSE])[,i]])
#         if(any(is.na(z[,i]))) z[,i] <- as.numeric((xz[,IND,drop=FALSE])[,i])
#       }
#     }
#     ## Don't assign names when factor.to.numeric is TRUE (otherwise
#     ## NAs populate matrix)
#     if(!factor.to.numeric) names(z) <- xznames[IND]
#     znames <- xznames[IND]
#     num.z <- ncol(z)
#   } else {
#     z <- NULL
#     znames <- NULL
#     num.z <- NULL
#   }
# 
#   return(list(x=x,
#               num.x=num.x,
#               xnames=xnames,
#               z=z,
#               num.z=num.z,
#               is.ordered.z=is.ordered.z,
#               znames=znames))
#   
# }

# Unit: microseconds
# expr     min       lq      mean   median       uq     max neval cld
# original 157.440 170.2730 185.73841 174.3730 183.3315 824.469   200  a
# optimized  71.832  78.6585  83.12197  81.5285  86.9815 110.372   200   b
# 
# -------- Benchmark Summary --------
# Median:  2.14x  ( 46.8% as fast)  - optimized * is fastest by median.
# Mean:    2.23x  ( 44.8% as fast)  - optimized * is fastest by mean.

# ---- Optimized (identical) ----
splitFrame <- function(xz, factor.to.numeric=FALSE) {
  if (missing(xz)) stop(" you must provide xz data")
  if (!is.data.frame(xz)) stop(" xz must be a data frame")
  xznames <- names(xz)
  IND <- vapply(xz, is.factor, logical(1))
  x <- xz[, !IND, drop = FALSE]
  num.x <- ncol(x)
  if (num.x == 0) stop(" can't fit spline surfaces with no continuous predictors")
  xnames <- xznames[!IND]
  is.ordered.z <- NULL
  
  if (any(IND)) {
    zdf <- xz[, IND, drop = FALSE]
    
    # IMPORTANT: match original exactly (no names)
    is.ordered.z <- unname(vapply(zdf, is.ordered, logical(1)))
    
    if (!factor.to.numeric) {
      z <- data.frame(zdf)
      names(z) <- xznames[IND]
    } else {
      # Plain numeric matrix, no dimnames (matches original loop)
      z <- sapply(zdf, function(col) {
        suppressWarnings(val <- as.numeric(levels(col))[col])
        if (any(is.na(val))) as.numeric(col) else val
      }, simplify = TRUE)
      
      if (is.vector(z)) z <- matrix(z, ncol = 1)
      storage.mode(z) <- "double"
      dimnames(z) <- NULL
    }
    
    znames <- xznames[IND]
    num.z <- if (!is.null(z)) NCOL(z) else 0L
  } else {
    z <- NULL
    znames <- NULL
    num.z <- NULL
  }
  
  list(
    x = x,
    num.x = num.x,
    xnames = xnames,
    z = z,
    num.z = num.z,
    is.ordered.z = is.ordered.z,
    znames = znames
  )
}


trim.quantiles = function(dat, trim){
  if (sign(trim) == sign(-1)){
    trim = abs(trim)
    tq = quantile(dat, probs = c(0.0, 0.0+trim, 1.0-trim,1.0))
    tq = c(2.0*tq[1]-tq[2], 2.0*tq[4]-tq[3])
  }
  else {
    tq = quantile(dat, probs = c(0.0+trim, 1.0-trim))
  }
  tq
}

uocquantile = function(x, prob) {
  if (is.ordered(x)){
    tq = unclass(table(x))
    tq = tq / sum(tq)
    j = which(sapply(1:length(tq), function(y){ sum(tq[1:y]) }) >= prob)[1]
    sort(unique(x))[j]
  } else if (is.factor(x)) {
    ## just returns mode
    tq = unclass(table(x))
    j = which(tq == max(tq))[1]
    sort(unique(x))[j]
  } else {
    quantile(x, probs = prob)
  }
}

## statistical functions

RSQfunc <- function(y,y.pred,weights=NULL) {
  if(!is.null(weights)) {
    y <- y*sqrt(weights)
    y.pred <- y.pred*sqrt(weights)
  }
  y.mean <- mean(y)
  return((sum((y-y.mean)*(y.pred-y.mean))^2)/(sum((y-y.mean)^2)*sum((y.pred-y.mean)^2)))
}

MSEfunc <- function(y,y.fit) {
  mean((y-y.fit)^2)
}

MAEfunc <- function(y,y.fit) {
  mean(abs(y-y.fit))
}

MAPEfunc <- function(y,y.fit) {
  jj = which(y != 0)
  
  mean(c(abs((y[jj]-y.fit[jj])/y[jj]), as.numeric(replicate(length(y)-length(jj),2))))
}

CORRfunc <- function(y,y.fit) {
  abs(corr(cbind(y,y.fit)))
}

SIGNfunc <- function(y,y.fit) {
  sum(sign(y) == sign(y.fit))/length(y)
}

# blank <- function(len){
#   sapply(len, function(nb){
#     paste(rep(' ', times = nb), collapse='')
#   })
# }

# Unit: milliseconds
# expr       min        lq      mean   median       uq      max neval cld
# original 320.34571 329.90273 339.40322 334.1152 346.4990 375.1420    50  a
# optimized  20.20201  20.26974  20.58103  20.4557  20.6688  23.3083    50   b
# 
# -------- Benchmark Summary --------
# Median:  16.33x  (  6.1% as fast)  - optimized * is fastest by median.
# Mean:    16.49x  (  6.1% as fast)  - optimized * is fastest by mean.

# ---- Optimized (identical) ----
blank <- function(len) {
  strrep(' ', len)
}

## regression quantile check function

# check.function <- function(u,tau=0.5) {
#   if(missing(u)) stop(" Error: u must be provided")
#   if(tau <= 0 | tau >= 1) stop(" Error: tau must lie in (0,1)")
#   return(u*(tau-ifelse(u<0,1,0)))
# }

# Unit: milliseconds
# expr       min        lq      mean    median        uq      max neval cld
# original 16.854977 17.637872 18.985206 17.985347 19.054730 30.00581   100  a
# optimized  2.483821  2.775454  3.372387  2.962947  3.430388 12.93263   100   b
# 
# -------- Benchmark Summary --------
# Median:  6.07x  ( 16.5% as fast)  - optimized * is fastest by median.
# Mean:    5.63x  ( 17.8% as fast)  - optimized * is fastest by mean.

# ---- Optimized (identical) ----
check.function <- function(u, tau = 0.5) {
  if (missing(u)) stop(" Error: u must be provided")
  if (tau <= 0 || tau >= 1) stop(" Error: tau must lie in (0,1)")
  u * (tau - (u < 0))
}


## Note - this is defined in cv.kernel.spline so if you modify there
## you must modify here also.

## Note - March 20 2012 - this is buggy - model$x is empty but
## hat(model$x) returns 1 so it passes. This is not used for
## cross-validation, rather only for summary/predict and potentially
## pruning, so for the moment we let it sit.

cv.rq <- function (model, tau = 0.5, weights = NULL) {
  return(mean(check.function(residuals(model),tau)/(1-hat(model$x))^(1/sqrt(tau*(1-tau)))))
}

## This function is based on functions in the limma package and
## corpcor package (is.positive.definite)... check the condition
## number of a matrix based on the ratio of max/min eigenvalue.  Note
## that the definition
## tol=max(dim(x))*max(sqrt(abs(e)))*.Machine$double.eps is exactly
## compatible with the conventions used in "Octave" or "Matlab".  Note
## that for weighted regression you simply use x*L which conducts
## row-wise multiplication (i.e. diag(L)%*%X not necessary). Note also
## that crossprod(X) is significantly faster than t(X)%*%X (matrix is
## symmetric so only use lower triangle).

is.fullrank <- function(x)
{
  e <- eigen(crossprod(as.matrix(x)), symmetric = TRUE, only.values = TRUE)$values
  e[1] > 0 && abs(e[length(e)]/e[1]) > max(dim(as.matrix(x)))*max(sqrt(abs(e)))*.Machine$double.eps
}

## Function that determines the dimension of the multivariate basis
## without precomputing it... the tensor is the mother that consumes
## ginormous amounts of memory, followed by the glp basis.

# dimBS <- function(basis="additive",kernel=TRUE,degree=NULL,segments=NULL,include=NULL,categories=NULL) {
# 
#   ## This function computes the dimension of the glp basis without the
#   ## memory overhead associated with computing the glp basis itself
#   ## (thanks to Zhenghua Nie)
# 
#   two.dimen<- function(d1,d2,nd1,pd12){
#     if(d2 ==1) {
#       ret <- list()
#       ret$d12 <- pd12
#       ret$nd1 <- nd1
#       return(ret)
#     }
#     d12 <- d2
#     if(d1-d2>0){
#       for(i in 1:(d1-d2)){
#         d12 <- d12+d2*nd1[i]
#       }}
#     if(d2>1){
#       for(i in 2:d2){
#         d12 <- d12 + (i*nd1[d1-i+1])
#       }
#     }
#     d12 <- d12 + nd1[d1]   ## The maximum number
#     
#     nd2 <- nd1  ## Calculate nd2
#     if(d1>1){
#       for(j in 1:(d1-1)) {
#         nd2[j] <- 0
#         for(i in j:max(0,j-d2+1)) {
#           if(i > 0) {
#             nd2[j] <- nd2[j] + nd1[i]                  
#           }
#           else {
#             nd2[j] <- nd2[j] + 1  ## nd1[0] always 1
#           }
#         }
#       }
#     }
#     if(d2>1) {
#       nd2[d1] <- nd1[d1]
#       for(i in (d1-d2+1):(d1-1)) nd2[d1] <- nd2[d1]+nd1[i]
#     }
#     else {
#       nd2[d1] <- nd1[d1]
#     }
#     ret <- list()
#     ret$d12 <- d12
#     ret$nd1 <- nd2 
#     
#     return(ret)
#   }
#   
#   ## Some basic error checking
#  
#   if(basis!="additive" & basis!="glp" & basis!="tensor") stop(" Error: basis must be either additive, glp, or tensor")
# 
#   if(!kernel)
#     if(is.null(include) | is.null(categories)) stop(" Error: you must provide include and categories vectors")    
# 
#   
#   K <- cbind(degree,segments)
# 
#   ncol.bs <- 0
# 
#   if(kernel) {
#     if(basis=="additive") {
#       if(any(K[,1] > 0))
#         ncol.bs <- sum(rowSums(K[K[,1]!=0,,drop=FALSE])-1)
#     }
#     if(basis=="glp") {
#       dimen <- rowSums(K[K[,1]!=0,,drop=FALSE])-1
#       dimen <- dimen[dimen>0] ## Delete elements which are equal to 0.
#       dimen <- sort(dimen,decreasing=TRUE) ## Sort the array to save memory when doing the computation.
#       k <-length(dimen)
#       if(k==0) {
#         ncol.bs <- 0
#       } else {
#         nd1 <- rep(1,dimen[1])   ## At the beginning,  we have one for [1, 2, 3, ..., dimen[1]]
#         nd1[dimen[1]] <- 0       ## nd1 represents the frequency for every element of [1, 2, 3, ..., dimen[1]]
#         ncol.bs <- dimen[1]
#         if(k>1) {
#           for(i in 2:k) {
#             dim.rt <- two.dimen(dimen[1],dimen[i],nd1,ncol.bs)
#             nd1 <- dim.rt$nd1
#             ncol.bs <- dim.rt$d12
#           }
#           ncol.bs <- dim.rt$d12+k-1
#         }
#       }
#     }
#     if(basis=="tensor") {
#       if(any(K[,1] > 0))
#         ncol.bs <- prod(rowSums(K[K[,1]!=0,,drop=FALSE]))
#     }
# 	} else {
# 			if(basis=="additive") {
# 					if(any(K[,1] > 0)) 
# 							ncol.bs <- sum(c(rowSums(K[K[,1]!=0,,drop=FALSE])-1,include*categories-1))
# 			}
# 			if(basis=="glp") {
# 					dimen <- c(rowSums(K[K[,1]!=0,,drop=FALSE])-1,include*categories-1)
# 					dimen <- dimen[dimen>0] ## Delete elements which are eqaul to 0.
# 					dimen <- sort(dimen,decreasing=TRUE) ## Sort the array to save memory when doing the computation.
# 					k <-length(dimen)
# 					if(k==0) {
# 							ncol.bs <- 0
# 					} else {
# 							nd1 <- rep(1,dimen[1])   ## At the beginning,  we have one for [1, 2, 3, ..., dimen[1]]
# 							nd1[dimen[1]] <- 0       ## nd1 represents the frequency for every element of [1, 2, 3, ..., dimen[1]]
# 							ncol.bs <- dimen[1]
# 							if(k>1) {
# 									for(i in 2:k) {
# 											dim.rt <- two.dimen(dimen[1],dimen[i],nd1,ncol.bs)
# 											nd1 <- dim.rt$nd1
# 											ncol.bs <- dim.rt$d12
# 									}
# 									ncol.bs <- dim.rt$d12+k-1
# 							}
# 					}
# 			}
# 			if(basis=="tensor") {
# 					if(any(K[,1] > 0)) 
# 							ncol.bs <- prod(c(rowSums(K[K[,1]!=0,,drop=FALSE]),(include*categories-1)))
# 			}
# 	}
# 
# 	return(ncol.bs)
# 
# }

# expr    min     lq     mean median     uq     max neval cld
# original 40.221 40.918 43.25090 41.369 43.870 120.663   200   a
# optimized 37.023 37.556 41.43398 38.089 41.164 377.938   200   a
# 
# -------- Benchmark Summary --------
# Median:  1.09x  ( 92.1% as fast)  - optimized * is fastest by median.
# Mean:    1.04x  ( 95.8% as fast)  - optimized * is fastest by mean.

# ---- Optimized (identical) ----
dimBS <- function(basis="additive", kernel=TRUE, degree=NULL, segments=NULL, include=NULL, categories=NULL) {
  two.dimen <- function(d1, d2, nd1, pd12) {
    if (d2 == 1) return(list(d12 = pd12, nd1 = nd1))
    d12 <- d2
    if (d1 > d2) {
      for (i in seq_len(d1 - d2)) d12 <- d12 + d2 * nd1[i]
    }
    if (d2 > 1) {
      for (i in 2:d2) d12 <- d12 + i * nd1[d1 - i + 1]
    }
    d12 <- d12 + nd1[d1]
    nd2 <- nd1
    if (d1 > 1) {
      for (j in 1:(d1 - 1)) {
        s <- 0
        for (i in j:max(0, j - d2 + 1)) s <- s + if (i > 0) nd1[i] else 1
        nd2[j] <- s
      }
    }
    if (d2 > 1) {
      nd2[d1] <- nd1[d1]
      for (i in (d1 - d2 + 1):(d1 - 1)) nd2[d1] <- nd2[d1] + nd1[i]
    } else {
      nd2[d1] <- nd1[d1]
    }
    list(d12 = d12, nd1 = nd2)
  }
  if (!basis %in% c('additive','glp','tensor')) stop(' Error: basis must be either additive, glp, or tensor')
  if (!kernel && (is.null(include) || is.null(categories))) stop(' Error: you must provide include and categories vectors')
  K <- cbind(degree, segments)
  ncol.bs <- 0
  if (kernel) {
    if (basis == 'additive') {
      if (any(K[,1] > 0)) ncol.bs <- sum(rowSums(K[K[,1] != 0, , drop=FALSE]) - 1)
    } else if (basis == 'glp') {
      rs <- rowSums(K[K[,1] != 0, , drop=FALSE]) - 1
      dimen <- sort(rs[rs > 0], decreasing = TRUE)
      k <- length(dimen)
      if (k == 0) {
        ncol.bs <- 0
      } else {
        nd1 <- rep(1, dimen[1]); nd1[dimen[1]] <- 0
        ncol.bs <- dimen[1]
        if (k > 1) {
          for (i in 2:k) {
            dim.rt <- two.dimen(dimen[1], dimen[i], nd1, ncol.bs)
            nd1 <- dim.rt$nd1
            ncol.bs <- dim.rt$d12
          }
          ncol.bs <- dim.rt$d12 + k - 1
        }
      }
    } else if (basis == 'tensor') {
      if (any(K[,1] > 0)) ncol.bs <- prod(rowSums(K[K[,1] != 0, , drop=FALSE]))
    }
  } else {
    if (basis == 'additive') {
      if (any(K[,1] > 0)) ncol.bs <- sum(c(rowSums(K[K[,1] != 0, , drop=FALSE]) - 1, include * categories - 1))
    } else if (basis == 'glp') {
      rs <- c(rowSums(K[K[,1] != 0, , drop=FALSE]) - 1, include * categories - 1)
      dimen <- sort(rs[rs > 0], decreasing = TRUE)
      k <- length(dimen)
      if (k == 0) {
        ncol.bs <- 0
      } else {
        nd1 <- rep(1, dimen[1]); nd1[dimen[1]] <- 0
        ncol.bs <- dimen[1]
        if (k > 1) {
          for (i in 2:k) {
            dim.rt <- two.dimen(dimen[1], dimen[i], nd1, ncol.bs)
            nd1 <- dim.rt$nd1
            ncol.bs <- dim.rt$d12
          }
          ncol.bs <- dim.rt$d12 + k - 1
        }
      }
    } else if (basis == 'tensor') {
      if (any(K[,1] > 0)) ncol.bs <- prod(c(rowSums(K[K[,1] != 0, , drop=FALSE]), (include * categories - 1)))
    }
  }
  ncol.bs
}
