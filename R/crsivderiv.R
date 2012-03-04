## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instruments
## xeval: optional evaluation data for the exogenous predictors

## ... optional arguments for crs()

## This function returns a list with the following elements:

## phi: the IV estimator of phi(z) corresponding to the estimated
## deriviative phi(z)

## This function will compute the integral using the trapezoidal rule
## and the cumsum function as we need to compute this in a
## computationally efficient manner.

integrate.trapezoidal <- function(x,y) {
  n <- length(x)
  rank.x <- rank(x)
  order.x <- order(x)
  y <- y[order.x]
  x <- x[order.x]
  int.vec <- numeric(length(x))
  int.vec[1] <- 0
  int.vec[2:n] <- cumsum((x[2:n] - x[2:n-1]) * (y[2:n] + y[2:n-1]) / 2)
  return(int.vec[rank.x])
}

## Avoid division by zero - for numerical integration etc.

NZD <- function(a) {
  sapply(1:NROW(a), function(i) {if(a[i] < 0) min(-.Machine$double.xmin,a[i]) else max(.Machine$double.xmin,a[i])})
}

crsivderiv <- function(y,
                       z,
                       w,
                       x=NULL,
                       zeval=NULL,
                       weval=NULL,
                       xeval=NULL,
                       iterate.max=100,
                       iterate.tol=1.0e-04,
                       constant=0.5,
                       start.phi.zero=FALSE,
                       stop.on.increase=TRUE,
                       smooth.residuals=TRUE,
                       opts=list("MAX_BB_EVAL"=10000,
                         "EPSILON"=.Machine$double.eps,
                         "INITIAL_MESH_SIZE"="r1.0e-01",
                         "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                         "MIN_POLL_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                         "DISPLAY_DEGREE"=0),
                       ...) {
  
  console <- newLineConsole()

  ## Basic error checking

  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")
  if(!is.logical(smooth.residuals)) stop("smooth.residuals must be logical (TRUE/FALSE)")  

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(!is.null(x) && NROW(y) != NROW(x)) stop("y and x have differing numbers of rows")
  if(iterate.max < 2) stop("iterate.max must be at least 2")

  ## Cast as data frames

  w <- data.frame(w)
  z <- data.frame(z)
  if(!is.null(x)) x <- data.frame(x)

  ## Check for evaluation data

  if(is.null(zeval)) zeval <- z
  if(is.null(weval)) weval <- w
  if(!is.null(x) && is.null(xeval)) xeval <- x

  ## Set up formulas for multivariate w, z, and x if provided

  wnames <- names(w)
  znames <- names(z)
  names(weval) <- wnames  
  names(zeval) <- znames  

  ## If there exist exogenous regressors X, append these to the
  ## formulas involving Z (can be manually added to W by the user if
  ## desired)

  if(!is.null(x)) {
    xnames <- names(x)
    names(xeval) <- xnames    
  }

  ## Now create evaluation data

  if(is.null(x)) {
    traindata <- data.frame(y,z,w)
    evaldata <- data.frame(zeval,weval)
  } else {
    traindata <- data.frame(y,z,w,x)    
    evaldata <- data.frame(zeval,weval,xeval)
  }

  ## Formulae for derivative estimation

  formula.muw <- as.formula(paste("mu ~ ", paste(wnames, collapse= "+")))
  formula.yw <- as.formula(paste("y ~ ", paste(wnames, collapse= "+")))
  formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+")))
  formula.phiw <- as.formula(paste("phi ~ ", paste(wnames, collapse= "+")))  

  ## Landweber-Fridman

  ## We begin the iteration computing phi.prime.0
    
  console <- printClear(console)
  console <- printPop(console)
  if(is.null(x)) {
    console <- printPush(paste("Computing optimal smoothing for f(z) and S(z) for iteration 1...",sep=""),console)
  } else {
    console <- printPush(paste("Computing optimal smoothing  f(z) and S(z) for iteration 1...",sep=""),console)
  }

  ## Note - here I am only treating the univariate case, so let's
  ## throw a stop with warning for now...

  if(NCOL(z) > 1) stop(" This version supports univariate z only")
  
  ## For all results we need the density function for Z and the
  ## survivor function for Z (1-CDF of Z)
  
  require(np)
  
  cat(paste("\rIteration ", 1, " of at most ", iterate.max,sep=""))

  ## Let's compute the bandwidth object for the unconditional
  ## density for the moment. Use the normal-reference rule for speed
  ## considerations.
  
  bw <- npudensbw(dat=z,bwmethod="normal-reference")
  model.fz <- npudens(tdat=z,bws=bw)
  f.z <- predict(model.fz,newdata=evaldata)
  model.Sz <- npudist(tdat=z,bws=bw)
  S.z <- 1-predict(model.Sz,newdata=evaldata)

  if(!start.phi.zero) {

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing for E(y|z) for iteration 1...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing  for E(y|z,x) for iteration 1...",sep=""),console)
    }
    
    model.E.y.z <- crs(formula.yz,
                       opts=opts,
                       data=traindata,
                       deriv=1,
                       ...)

    E.y.z <- predict(model.E.y.z,newdata=evaldata)
    
    phi.prime <- attr(E.y.z,"deriv.mat")[,1]
    
    ## Step 1 - begin iteration - for this we require \varphi_0. To
    ## compute \varphi_{0,i}, we require \mu_{0,i}. For j=0 (first
    ## term in the series), \mu_{0,i} is Y_i.

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing for E(y|w) (stopping rule) for iteration 1...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing  for E(y|w) (stopping rule) for iteration 1...",sep=""),console)
    }

    ## For stopping rule...

    model.E.y.w <- crs(formula.yw,
                      opts=opts,
                      data=traindata,
                      ...)

    E.y.w <- predict(model.E.y.w,newdata=evaldata)

  } else {
    
    ## Step 1 - begin iteration - for this we require \varphi_0. To
    ## compute \varphi_{0,i}, we require \mu_{0,i}. For j=0 (first
    ## term in the series), \mu_{0,i} is Y_i.
    
    mu <- y
    
    ## We also require the mean of \miu_{0,i}
    
    mean.mu <- mean(mu)
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing for E(y|w) (stopping rule) for iteration 1...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing  for E(y|w) (stopping rule) for iteration 1...",sep=""),console)
    }

    ## Next, we regress require \mu_{0,i} W

    model.E.y.w <- crs(formula.yw,
                      opts=opts,
                      data=traindata,
                      ...)

    E.y.w <- predict(model.E.y.w,newdata=evaldata)

    ## We require the mean of the fitted values
    
    mean.predicted.E.mu.w <- mean(E.y.w)
    
    ## We need the mean of the fitted values for this (we readily
    ## compute the CDF not the survivor, so anything that is weighted
    ## by the survivor kernel can be expressed as the mean of that
    ## being weighted minus the weighting using the CDF kernel).
    ## Next, we need the weighted sum of the survivor kernel where the
    ## weights are E[\mu_{0,i}|W]. We can write this as the mean of the
    ## \mu_{0,i} minus the weighted sum using the CDF kernel, i.e. if
    ## K is a CDF kernel, then n^{-1}\sum_j \bar K() \mu_{0,i} =
    ## n^{-1}\sum_j (1- K()) \mu_{0,i} = n^{-1}\sum_j\mu_{0,i}-
    ## n^{-1}\sum_j K() \mu_{0,i}
    
    ## Now we compute T^* applied to E.y.w, and this is phi.prime.0 for
    ## j=0.
    
    ## CDF weighted sum (but we need survivor weighted sum...)
    
    cdf.weighted.average <- npksum(txdat=z,
                                   exdat=zeval,
                                   tydat=as.matrix(E.y.w),
                                   operator="integral",
                                   bws=bw$bw)$ksum/length(y)
    
    survivor.weighted.average <- mean.predicted.E.mu.w - cdf.weighted.average

    phi.prime <- (survivor.weighted.average - S.z*mean.mu)/f.z

  }  
      
  norm.stop <- numeric()

  ## NOTE - this presumes univariate z case... in general this would
  ## be a continuous variable's index

  phi <- integrate.trapezoidal(z[,1],phi.prime)
  
  ## In the definition of phi we have the integral minus the mean of
  ## the integral with respect to z, so subtract the mean here
  
  phi <- phi - mean(phi) + mean(y)
  
  ## For the stopping rule, we require E.phi.w
  
  model.E.phi.w <- crs(formula.phiw,
                       opts=opts,
                       data=traindata,
                       ...)
  
  E.phi.w <- predict(model.E.phi.w,newdata=evaldata)

  norm.stop[1] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)
  
  ## Now we compute mu.0 (a residual of sorts)
  
  mu <- y - phi
  
  ## Now we repeat this entire process using mu = y = phi.0 rather than y
  
  mean.mu <- mean(mu)
  
  ## Next, we regress require \mu_{0,i} W
  
  if(smooth.residuals) {

    ## Smooth residuals

    model.E.mu.w <- crs(formula.muw,
                        opts=opts,data=traindata,
                        cv="none",
                        degree=model.E.phi.w$degree,
                        segments=model.E.phi.w$segments,
                        ...)
    
    ## We require the fitted values...
    
    predicted.model.E.mu.w <- predict(model.E.mu.w,newdata=evaldata)
    
    ## We again require the mean of the fitted values
    
    mean.predicted.model.E.mu.w <- mean(predicted.model.E.mu.w)
    
  } else {

    model.E.phi.w <- crs(formula.phiw,
                        opts=opts,data=traindata,
                        cv="none",
                        degree=model.E.phi.w$degree,
                        segments=model.E.phi.w$segments,
                        ...)

    ## We require the fitted values...
    
    predicted.model.E.phi.w <- E.y.w - predict(model.E.phi.w,newdata=evaldata)
    
    ## We again require the mean of the fitted values
    
    mean.predicted.model.E.phi.w <- mean(E.y.w) - mean(predicted.model.E.phi.w)
    
  }
  
  ## Now we compute T^* applied to mu
  
  cdf.weighted.average <- npksum(txdat=z,
                                 exdat=zeval,
                                 tydat=as.matrix(predicted.model.E.mu.w),
                                 operator="integral",
                                 bws=bw$bw)$ksum/nrow(traindata)
  
  survivor.weighted.average <- mean.predicted.model.E.mu.w - cdf.weighted.average
  
  T.star.mu <- (survivor.weighted.average-S.z*mean.mu)/f.z
    
  ## Now we update phi.prime.0, this provides phi.prime.1, and now
  ## we can iterate until convergence... note we replace phi.prime.0
  ## with phi.prime.1 (i.e. overwrite phi.prime)
  
  phi.prime <- phi.prime + constant*T.star.mu
  
  ## This we iterate...
  
  for(j in 2:iterate.max) {

    ## Save previous run in case stop norm increases
    
    phi.j.m.1 <- phi
    phi.prime.j.m.1 <- phi.prime
    
    cat(paste("\rIteration ", j, " of at most ", iterate.max,sep=""))
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration ", j,"...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing and phi(z,x) for iteration ", j,"...",sep=""),console)
    }
    
    ## NOTE - this presumes univariate z case... in general this would
    ## be a continuous variable's index
    
    phi <- integrate.trapezoidal(z[,1],phi.prime)
    
    ## In the definition of phi we have the integral minus the mean of
    ## the integral with respect to z, so subtract the mean here
    
    phi <- phi - mean(phi) + mean(y)
    
    ## For the stopping rule, we require E.phi.w
    
    model.E.phi.w <- crs(formula.phiw,
                         opts=opts,
                         data=traindata,
                         ...)
    
    E.phi.w <- predict(model.E.phi.w,newdata=evaldata)
    norm.stop[j] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)
    
    ## Now we compute mu.0 (a residual of sorts)
    
    mu <- y - phi
    
    ## Now we repeat this entire process using mu = y = phi.0 rather than y
    
    mean.mu <- mean(mu)
    
    ## Next, we regress require \mu_{0,i} W

    if(smooth.residuals) {
      
      ## Smooth residuals
      
      model.E.mu.w <- crs(formula.muw,
                          opts=opts,data=traindata,
                          cv="none",
                          degree=model.E.phi.w$degree,
                          segments=model.E.phi.w$segments,
                          ...)
      
      ## We require the fitted values...
      
      predicted.model.E.mu.w <- predict(model.E.mu.w,newdata=evaldata)
      
      ## We again require the mean of the fitted values
      
      mean.predicted.model.E.mu.w <- mean(predicted.model.E.mu.w)
      
    } else {
      
      model.E.phi.w <- crs(formula.phiw,
                           opts=opts,data=traindata,
                           cv="none",
                           degree=model.E.phi.w$degree,
                           segments=model.E.phi.w$segments,
                           ...)
      
      ## We require the fitted values...
      
      predicted.model.E.my.w <- E.y.w - predict(model.E.phi.w,newdata=evaldata)
      
      ## We again require the mean of the fitted values
      
      mean.predicted.model.E.phi.w <- mean(E.y.w) - mean(predicted.model.E.phi.w)
      
    }
    
    ## Now we compute T^* applied to mu
    
    cdf.weighted.average <- npksum(txdat=z,
                                   exdat=zeval,
                                   tydat=as.matrix(predicted.model.E.mu.w),
                                   operator="integral",
                                   bws=bw$bw)$ksum/nrow(traindata)
    
    survivor.weighted.average <- mean.predicted.model.E.mu.w - cdf.weighted.average
    
    T.star.mu <- (survivor.weighted.average-S.z*mean.mu)/f.z
    
    ## Now we update, this provides phi.prime.1, and now we can iterate until convergence...
      
    phi.prime <- phi.prime + constant*T.star.mu

    ## If stopping rule criterion increases or we are below stopping
    ## tolerance then break
    
    if(norm.stop[j] < iterate.tol) break()
    if(stop.on.increase && norm.stop[j] > norm.stop[j-1]) {
      phi <- phi.j.m.1 
      phi.prime <- phi.prime.j.m.1
      break()
    }
    
  }

  console <- printClear(console)
  console <- printPop(console)
  
  if(j == iterate.max) warning(" iterate.max reached: increase iterate.max or inspect norm.stop vector")
  
  return(list(phi=phi,phi.prime=phi.prime,num.iterations=j,norm.stop=norm.stop))
  
}
