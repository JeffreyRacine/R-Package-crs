## This function tests for the maximum well-conditioned spline degree.

## Note that increasing the number of breaks, other things equal,
## results in a better-conditioned matrix. Hence we ignore nbreak and
## set it to its minimum (2)

check.max.spline.degree <- function(xdat=NULL,degree=NULL,issue.warning=FALSE) {

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
        if(rcond(t(X)%*%X)<.Machine$double.eps) {
          for(j in 1:degree[i]) {
            d[i] <- j
            X <- gsl.bs(xdat[,numeric.index[i]],degree=d[i],nbreak=2)
            if(rcond(t(X)%*%X)<.Machine$double.eps) {
              d[i] <- j-1
              break()
            }
          }
        }
        if(d[i] < degree[i]) {
          if(issue.warning) warning(paste("\r Predictor ",i," B-spline basis is ill-conditioned beyond degree ",d[i],": see note in ?npglpreg",sep=""),immediate.=TRUE)
          ill.conditioned <- TRUE
        }
      }
    }

  }

  attr(ill.conditioned, "degree.max.vec") <- d
  return(ill.conditioned)

}

## Utility function for dimension of par(mfrow=c(,)) for multiple
## plots on the same device

dim.plot = function(x) {
  a1 = round(sqrt(4.0/3.0*x))
  a2 = ceiling(x/a1)
  c(a1,a2)
}

succeedWithResponse <- function(tt, frame){
  !any(class(try(eval(expr = attr(tt, "variables"),
                      envir = frame, encl = NULL), silen = TRUE)) == "try-error")
}

## Utility function to divide explanatory variables into
## factors/numeric, strip off names etc.

splitFrame <- function(xz, factor.to.numeric=FALSE) {
  
  if(missing(xz)) stop(" you must provide xz data")
  if(!is.data.frame(xz)) stop(" xz must be a data frame")

  xznames <- names(xz)
  
  IND <- logical()

  for(i in 1:NCOL(xz)) IND[i] <- is.factor(xz[,i])

  x <- xz[,!IND,drop=FALSE]
  num.x <- ncol(x)

  ## We require at least one continuous predictor to conduct spline
  ## smoothing, but there may/may not be factors.

  if(num.x == 0) stop(" can't fit spline surfaces with no continuous predictors")

  xnames <- xznames[!IND]

  is.ordered.z <- NULL
  
  if(any(IND)) {
    is.ordered.z <- logical()
    for(i in 1:NCOL(xz[,IND,drop=FALSE])) is.ordered.z[i] <- is.ordered((xz[,IND,drop=FALSE])[,i])
    if(!factor.to.numeric) {
      z <- data.frame(xz[,IND,drop=FALSE])
    } else {
      ## If factor.to.numeric crudely convert factors to numeric.
      z <- matrix(NA,NROW(xz),NCOL(xz[,IND,drop=FALSE]))
      ## To “revert” a factor f to its original numeric values,
      ## as.numeric(levels(f))[f] is recommended. Problem is that for
      ## character strings it produces a warning message. No idea how
      ## to test for this so dropping for the moment. Will affect
      ## ordered types.
      for(i in 1:NCOL(xz[,IND,drop=FALSE])) {
        suppressWarnings(z[,i] <- as.numeric(levels((xz[,IND,drop=FALSE])[,i]))[(xz[,IND,drop=FALSE])[,i]])
        if(any(is.na(z[,i]))) z[,i] <- as.numeric((xz[,IND,drop=FALSE])[,i])
      }
    }
    ## Don't assign names when factor.to.numeric is TRUE (otherwise
    ## NAs populate matrix)
    if(!factor.to.numeric) names(z) <- xznames[IND]
    znames <- xznames[IND]
    num.z <- ncol(z)
  } else {
    z <- NULL
    znames <- NULL
    num.z <- NULL
  }

  return(list(x=x,
              num.x=num.x,
              xnames=xnames,
              z=z,
              num.z=num.z,
              is.ordered.z=is.ordered.z,
              znames=znames))
  
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

RSQfunc <- function(y,y.pred) {
  y.mean <- mean(y)
  (sum((y-y.mean)*(y.pred-y.mean))^2)/(sum((y-y.mean)^2)*sum((y.pred-y.mean)^2))
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

blank <- function(len){
  sapply(len, function(nb){
    paste(rep(' ', times = nb), collapse='')
  })
}

## regression quantile check function

check.function <- function(u,tau=0.5) {
  if(missing(u)) stop(" Error: u must be provided")
  if(tau <= 0 | tau >= 1) stop(" Error: tau must lie in (0,1)")
  return(u*(tau-ifelse(u<0,1,0)))
}

## Note - this is defined in cv.kernel.spline so if you modify there
## you must modify here also

cv.rq <- function (model, tau = 0.5) {
  return(mean(check.function(residuals(model),tau)/(1-hat(model$x))^(1/sqrt(tau*(1-tau)))))
}

## From the limma package... check the condition number of a matrix
## (ratio of max/min eigenvalue) using .Machine$double.eps rather than
## their 1e-13 constant. Note that for weighted regression you simply
## use x*L which conducts row-wise multiplication (i.e. diag(L)%*%X
## not necessary)

is.fullrank <- function (x) 
{
    x <- as.matrix(x)
    e <- eigen(crossprod(x), symmetric = TRUE, only.values = TRUE)$values
    e[1] > 0 && abs(e[length(e)]/e[1]) > .Machine$double.eps
}
