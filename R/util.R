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

  if(any(IND)) {
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
      for(i in 1:NCOL(xz[,IND,drop=FALSE])) z[,i] <- as.numeric((xz[,IND,drop=FALSE])[,i])
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
              znames=znames))
  
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
