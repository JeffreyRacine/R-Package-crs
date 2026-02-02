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
## derivative phihat(z)
## phi.prime: the IV derivative estimator
## phi.mat: the matrix with colums phi_1, phi_2 etc. over all iterations
## phi.prime.mat: the matrix with colums phi'_1, phi'_2 etc. over all iterations
## num.iterations: number of iterations taken by Landweber-Fridman
## norm.stop: the stopping rule for each Landweber-Fridman iteration
## norm.value: the norm not multiplied by the number of iterations
## convergence: a character string indicating whether/why iteration terminated

crsivderiv <- function(y,
                       z,
                       w,
                       x=NULL,
                       zeval=NULL,
                       weval=NULL,
                       xeval=NULL,
                       constant=0.5,
                       display.nomad.progress=TRUE,
                       display.warnings=TRUE,
                       iterate.diff.tol=1.0e-08,
                       iterate.max=1000,
                       opts=list("MAX_BB_EVAL"=10000,
                                 "EPSILON"=.Machine$double.eps,
                                 "INITIAL_MESH_SIZE"="r1.0e-01",
                                 "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                                 "MIN_POLL_SIZE"=paste("r",1,sep=""),
                                 "DISPLAY_DEGREE"=0),
                       penalize.iteration=TRUE,
                       smooth.residuals=TRUE,
                       start.from=c("Eyz","EEywz"),
                       starting.values=NULL,
                       stop.on.increase=TRUE,
                       ...) {
  
  ptm.start <- proc.time()
  crs.messages <- getOption("crs.messages")
  is.eval.train <- is.null(zeval) && is.null(weval) && is.null(xeval)
  
  dot.args <- list(...)
  nmulti.user <- dot.args$nmulti
  nmulti <- if(!is.null(nmulti.user)) nmulti.user else 5
  
  ## Strip arguments we might pass explicitly to avoid conflicts
  dots.no.nmulti <- dot.args
  dots.no.nmulti$nmulti <- NULL
  
  ## Create a version of dots for pre-loop calls that strips conflicting args
  dots.preloop <- dot.args
  dots.preloop$formula <- NULL
  dots.preloop$opts <- NULL
  dots.preloop$data <- NULL
  dots.preloop$display.nomad.progress <- NULL
  dots.preloop$display.warnings <- NULL
  
  ## Create a version of dots for the loop that excludes parameters we pass explicitly
  dots.loop <- dots.no.nmulti
  dots.loop$formula <- NULL
  dots.loop$opts <- NULL
  dots.loop$data <- NULL
  dots.loop$display.nomad.progress <- NULL
  dots.loop$display.warnings <- NULL
  dots.loop$degree <- NULL
  dots.loop$segments <- NULL
  dots.loop$lambda <- NULL
  dots.loop$include <- NULL
  dots.loop$nmulti <- NULL # already handled in dots.no.nmulti but for safety
  
  ## Determine nmulti for subsequent loop iterations
  nmulti.loop <- if(!is.null(nmulti.user)) nmulti.user else 1
  
  console <- newLineConsole()
  
  ## Basic error checking
  
  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")
  if(!is.logical(smooth.residuals)) stop("smooth.residuals must be logical (TRUE/FALSE)")
  start.from <- match.arg(start.from)
  
  if(constant <= 0 || constant >=1) stop("constant must lie in (0,1)")
  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(!is.null(x) && NROW(y) != NROW(x)) stop("y and x have differing numbers of rows")
  if(iterate.max < 2) stop("iterate.max must be at least 2")
  if(iterate.diff.tol < 0) stop("iterate.diff.tol must be non-negative")
  
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
  
  if(!is.null(starting.values) && (NROW(starting.values) != NROW(evaldata))) stop(paste("starting.values must be of length",NROW(evaldata)))
  
  ## Formulae for derivative estimation
  
  formula.muw <- as.formula(paste("mu ~ ", paste(wnames, collapse= "+")))
  formula.yw <- as.formula(paste("y ~ ", paste(wnames, collapse= "+")))
  formula.phiw <- as.formula(paste("phi ~ ", paste(wnames, collapse= "+")))
  if(is.null(x)) {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+")))
  } else {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
  }
  
  ## Landweber-Fridman
  
  ## We begin the iteration computing phi.prime.0
  
  console <- printClear(console)
  console <- printPop(console)
  if(is.null(x)) {
    if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing for f(z) and S(z) for iteration 1"," of at most ", iterate.max,"...",sep=""),console)
  } else {
    if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing  f(z) and S(z) for iteration 1"," of at most ", iterate.max,"...",sep=""),console)
  }
  
  ## Note - here I am only treating the univariate case, so let's
  ## throw a stop with warning for now...
  
  if(NCOL(z) > 1) stop(" This version supports univariate z only")
  
  ## For all results we need the density function for Z and the
  ## survivor function for Z (1-CDF of Z)
  
  #  require(np)
  
  ## Let's compute the bandwidth object for the unconditional density
  ## for the moment. Use the normal-reference rule for speed
  ## considerations (sensitivity analysis indicates this is not
  ## problematic).
  
  bw <- npudensbw(dat=z,
                  bwmethod="normal-reference",
                  ...)
  model.fz <- npudens(tdat=z,
                      bws=bw$bw,
                      ...)
  f.z <- if(is.eval.train) predict(model.fz) else predict(model.fz,newdata=evaldata)
  model.Sz <- npudist(tdat=z,
                      bws=bw$bw,
                      ...)
  S.z <- 1-(if(is.eval.train) predict(model.Sz) else predict(model.Sz,newdata=evaldata))
  
  if(is.null(starting.values)) {
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing for E(y|z) for iteration 1 of at most ", iterate.max,"...",sep=""),console)
    } else {
      if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing  for E(y|z,x) for iteration 1 of at most ", iterate.max,"...",sep=""),console)
    }
    
    if(start.from == "Eyz") {
      ## Start from E(Y|z)      
      if(crs.messages) options(crs.messages=FALSE)
      model.E.y.z <- crs(formula.yz,
                         opts=opts,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,
                         data=traindata,
                         deriv=1,
                         ...)
      if(crs.messages) options(crs.messages=TRUE)    
      phi.0 <- model.E.y.z
      if(is.eval.train) {
        E.y.z <- model.E.y.z$fitted.values
        phi.prime <- model.E.y.z$deriv.mat[,1]
      } else {
        E.y.z_tmp <- predict(model.E.y.z,newdata=evaldata)
        E.y.z <- as.numeric(E.y.z_tmp)
        phi.prime <- attr(E.y.z_tmp,"deriv.mat")[,1]
      }
    } else {
      ## Start from E(E(Y|w)|z)
      tmp.model <- crs(formula.yw,opts=opts,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,data=traindata,...)
      E.y.w <- fitted(tmp.model)
      model.E.E.y.w.z <- crs(formula.Eywz,opts=opts,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,data=traindata,deriv=1,...)
      phi.0 <- model.E.E.y.w.z
      if(is.eval.train) {
        E.E.y.w.z <- model.E.E.y.w.z$fitted.values
        phi.prime <- model.E.E.y.w.z$deriv.mat[,1]
      } else {
        E.E.y.w.z_tmp <- predict(model.E.E.y.w.z,newdata=evaldata,...)
        E.E.y.w.z <- as.numeric(E.E.y.w.z_tmp)
        phi.prime <- attr(E.E.y.w.z_tmp,"deriv.mat")[,1]
      }
    }
    
  } else {
    phi.prime <- starting.values
    if(crs.messages) options(crs.messages=FALSE)
    phi.0 <- crs(formula.yz,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
    if(crs.messages) options(crs.messages=TRUE)
  }
  
  ## Step 1 - begin iteration - for this we require \varphi_0. To
  ## compute \varphi_{0,i}, we require \mu_{0,i}. For j=0 (first
  ## term in the series), \mu_{0,i} is Y_i.
  
  console <- printClear(console)
  console <- printPop(console)
  if(is.null(x)) {
    if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing for E(y|w) (stopping rule) for iteration 1 of at most ", iterate.max,"...",sep=""),console)
  } else {
    if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing  for E(y|w) (stopping rule) for iteration 1 of at most ", iterate.max,"...",sep=""),console)
  }
  
  ## NOTE - this presumes univariate z case... in general this would
  ## be a continuous variable's index
  
  phi <- integrate.trapezoidal(z[,1],phi.prime)
  
  ## In the definition of phi we have the integral minus the mean of
  ## the integral with respect to z, so subtract the mean here
  
  phi <- phi - mean(phi) + mean(y)
  
  starting.values.phi <- phi
  starting.values.phi.prime <- phi.prime
  
  ## For stopping rule...
  
  if(crs.messages) options(crs.messages=FALSE)
  model.E.y.w <- crs(formula.yw,
                     opts=opts,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,
                     data=traindata,
                     ...)
  if(crs.messages) options(crs.messages=TRUE)    
  
  E.y.w <- predict(model.E.y.w,newdata=evaldata)
  
  norm.stop <- numeric()
  
  ## For the stopping rule, we require E.phi.w
  
  if(crs.messages) options(crs.messages=FALSE)
  traindata$phi <- phi
  model.E.phi.w <- crs(formula.phiw,
                       opts=opts,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,
                       data=traindata,
                       ...)
  if(crs.messages) options(crs.messages=TRUE)    
  
  E.phi.w <- if(is.eval.train) fitted(model.E.phi.w) else predict(model.E.phi.w,newdata=evaldata)
  
  ## Now we compute mu.0 (a residual of sorts)
  
  mu <- y - phi
  
  ## Now we repeat this entire process using mu = y - phi.0 rather
  ## than y
  
  mean.mu <- mean(mu)
  
  if(smooth.residuals) {
    
    ## Smooth residuals (smooth of (y-phi) on w)
    
    if(crs.messages) options(crs.messages=FALSE)
    traindata$mu <- mu
    
    model.E.mu.w <- do.call("crs", c(list(formula=formula.muw,
                                          opts=opts,
                                          display.nomad.progress=display.nomad.progress,
                                          display.warnings=display.warnings,
                                          data=traindata),
                                     dots.preloop))
    
    ## Capture initial parameters for warm start
    degree.muw <- model.E.mu.w$degree
    segments.muw <- model.E.mu.w$segments
    lambda.muw <- model.E.mu.w$lambda
    include.muw <- model.E.mu.w$include
    
    ## Initialize unused warm start parameters
    degree.phiw <- NULL
    segments.phiw <- NULL
    lambda.phiw <- NULL
    include.phiw <- NULL
    
    if(crs.messages) options(crs.messages=TRUE)    
    
    ## We require the fitted values...
    
    predicted.model.E.mu.w <- if(is.eval.train) fitted(model.E.mu.w) else predict(model.E.mu.w,newdata=evaldata)
    
    ## We again require the mean of the fitted values
    
    mean.predicted.model.E.mu.w <- mean(predicted.model.E.mu.w)
    
  } else {
    
    ## Not smoothing residuals (difference of E(Y|w) and smooth of phi
    ## on w)
    
    if(crs.messages) options(crs.messages=FALSE)
    traindata$phi <- phi
    
    model.E.phi.w <- do.call("crs", c(list(formula=formula.phiw,
                                           opts=opts,
                                           display.nomad.progress=display.nomad.progress,
                                           display.warnings=display.warnings,
                                           data=traindata),
                                     dots.preloop))
    
    ## Capture initial parameters for warm start
    degree.phiw <- model.E.phi.w$degree
    segments.phiw <- model.E.phi.w$segments
    lambda.phiw <- model.E.phi.w$lambda
    include.phiw <- model.E.phi.w$include
    
    ## Initialize unused warm start parameters
    degree.muw <- NULL
    segments.muw <- NULL
    lambda.muw <- NULL
    include.muw <- NULL
    
    if(crs.messages) options(crs.messages=TRUE)    
    
    ## We require the fitted values...
    
    predicted.model.E.mu.w <- E.y.w - (if(is.eval.train) fitted(model.E.phi.w) else predict(model.E.phi.w,newdata=evaldata))
    
    ## We again require the mean of the fitted values
    
    mean.predicted.model.E.mu.w <- mean(E.y.w) - mean(predicted.model.E.mu.w)
    
  }
  
  norm.stop[1] <- sum(predicted.model.E.mu.w^2)/sum(E.y.w^2)
  
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
  
  phi.prime.mat <- matrix(NA, nrow=length(phi.prime), ncol=iterate.max)
  phi.mat <- matrix(NA, nrow=length(phi), ncol=iterate.max)
  phi.prime.mat[,1] <- phi.prime
  phi.mat[,1] <- phi
  
  ## Create a version of dots for the loop that excludes parameters we pass explicitly
  ## to avoid "matched by multiple actual arguments" errors.
  dots.loop <- dots.no.nmulti
  dots.loop$degree <- NULL
  dots.loop$segments <- NULL
  dots.loop$lambda <- NULL
  dots.loop$include <- NULL
  
  ## This we iterate...
  
  for(j in 2:iterate.max) {
    
    ## Save previous run in case stop norm increases
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration ", j," of at most ", iterate.max,"...",sep=""),console)
    } else {
      if(display.nomad.progress) console <- printPush(paste("Computing optimal smoothing and phi(z,x) for iteration ", j," of at most ", iterate.max,"...",sep=""),console)
    }
    
    ## NOTE - this presumes univariate z case... in general this would
    ## be a continuous variable's index
    
    phi <- integrate.trapezoidal(z[,1],phi.prime)
    
    ## In the definition of phi we have the integral minus the mean of
    ## the integral with respect to z, so subtract the mean here
    
    phi <- phi - mean(phi) + mean(y)
    
    ## Now we compute mu.0 (a residual of sorts)
    
    mu <- y - phi
    
    ## Now we repeat this entire process using mu = y = phi.0 rather than y
    ## Next, we regress require \mu_{0,i} W
    
    if(smooth.residuals) {
      
      ## Smooth residuals (smooth of (y-phi) on w)
      
      if(crs.messages) options(crs.messages=FALSE)
      traindata$mu <- mu
      
      model.E.mu.w <- do.call("crs", c(list(formula=formula.muw,
                                              degree=degree.muw,
                                              segments=segments.muw,
                                              lambda=lambda.muw,
                                              include=include.muw,
                                              nmulti=nmulti.loop,
                                              opts=opts,
                                              display.nomad.progress=display.nomad.progress,
                                              display.warnings=display.warnings,
                                              data=traindata),
                                         dots.loop))
      
      degree.muw <- model.E.mu.w$degree
      segments.muw <- model.E.mu.w$segments
      lambda.muw <- model.E.mu.w$lambda
      include.muw <- model.E.mu.w$include
      
      ## We require the fitted values...
      
      predicted.model.E.mu.w <- if(is.eval.train) fitted(model.E.mu.w) else predict(model.E.mu.w,newdata=evaldata)
      if(crs.messages) options(crs.messages=TRUE)    
      
      ## We again require the mean of the fitted values
      
      mean.predicted.model.E.mu.w <- mean(predicted.model.E.mu.w)
      
    } else {
      
      ## Not smoothing residuals (difference of E(Y|w) and smooth of
      ## phi on w)
      
      if(crs.messages) options(crs.messages=FALSE)
      traindata$phi <- phi
      
      model.E.phi.w <- do.call("crs", c(list(formula=formula.phiw,
                                               degree=degree.phiw,
                                               segments=segments.phiw,
                                               lambda=lambda.phiw,
                                               include=include.phiw,
                                               nmulti=nmulti.loop,
                                               opts=opts,
                                               display.nomad.progress=display.nomad.progress,
                                               display.warnings=display.warnings,
                                               data=traindata),
                                          dots.loop))
      
      degree.phiw <- model.E.phi.w$degree
      segments.phiw <- model.E.phi.w$segments
      lambda.phiw <- model.E.phi.w$lambda
      include.phiw <- model.E.phi.w$include
      
      ## We require the fitted values...
      
      predicted.model.E.mu.w <- E.y.w - (if(is.eval.train) fitted(model.E.phi.w) else predict(model.E.phi.w,newdata=evaldata))
      if(crs.messages) options(crs.messages=TRUE)    
      
      ## We again require the mean of the fitted values
      
      mean.predicted.model.E.mu.w <- mean(E.y.w) - mean(predicted.model.E.mu.w)
      
    }
    
    norm.stop[j] <- ifelse(penalize.iteration,j*sum(predicted.model.E.mu.w^2)/sum(E.y.w^2),sum(predicted.model.E.mu.w^2)/sum(E.y.w^2))
    
    ## Now we compute T^* applied to mu
    
    cdf.weighted.average <- npksum(txdat=z,
                                   exdat=zeval,
                                   tydat=as.matrix(predicted.model.E.mu.w),
                                   operator="integral",
                                   bws=bw$bw)$ksum/nrow(traindata)
    
    survivor.weighted.average <- mean.predicted.model.E.mu.w - cdf.weighted.average
    
    T.star.mu <- (survivor.weighted.average-S.z*mean.predicted.model.E.mu.w)/f.z
    
    ## Now we update, this provides phi.prime.1, and now we can iterate until convergence...
    
    phi.prime <- phi.prime + constant*T.star.mu
    phi.prime.mat[,j] <- phi.prime
    phi.mat[,j] <- phi
    
    ## The number of iterations in LF is asymptotically equivalent to
    ## 1/alpha (where alpha is the regularization parameter in
    ## Tikhonov).  Plus the criterion function we use is increasing
    ## for very small number of iterations. So we need a threshold
    ## after which we can pretty much confidently say that the
    ## stopping criterion is decreasing.  In Darolles et al. (2011)
    ## \alpha ~ O(N^(-1/(min(beta,2)+2)), where beta is the so called
    ## qualification of your regularization method. Take the worst
    ## case in which beta = 0 and then the number of iterations is ~
    ## N^0.5. Note that derivative estimation seems to require more
    ## iterations hence the heuristic sqrt(N)
    
    if(j > round(sqrt(nrow(traindata)))  && !is.monotone.increasing(norm.stop)) {
      ## If stopping rule criterion increases or we are below stopping
      ## tolerance then break
      
      if(stop.on.increase && norm.stop[j] > norm.stop[j-1]) {
        convergence <- "STOP_ON_INCREASE"
        break()
      }
      if(abs(norm.stop[j-1]-norm.stop[j]) < iterate.diff.tol) {
        convergence <- "ITERATE_DIFF_TOL"
        break()
      }
      
    }
    
    convergence <- "ITERATE_MAX"
    
  }
  
  phi.mat <- phi.mat[, 1:length(norm.stop), drop = FALSE]
  phi.prime.mat <- phi.prime.mat[, 1:length(norm.stop), drop = FALSE]
  
  ## Extract minimum, and check for monotone increasing function and
  ## issue warning in that case. Otherwise allow for an increasing
  ## then decreasing (and potentially increasing thereafter) portion
  ## of the stopping function, ignore the initial increasing portion,
  ## and take the min from where the initial inflection point occurs
  ## to the length of norm.stop
  
  norm.value <- norm.stop/(1:length(norm.stop))
  
  if(which.min(norm.stop) == 1 && is.monotone.increasing(norm.stop)) {
    if(display.warnings) warning("Stopping rule increases monotonically (consult model$norm.stop):\nThis could be the result of an inspired initial value (unlikely)\nNote: we suggest manually choosing phi.0 and restarting (e.g. instead set `starting.values' to E[E(Y|w)|z])")
    convergence <- "FAILURE_MONOTONE_INCREASING"
    #    phi <- starting.values.phi
    #    phi.prime <- starting.values.phi.prime
    j <- 1
    while(j < length(norm.value) && norm.value[j+1] > norm.value[j]) j <- j + 1
    j <- j-1 + which.min(norm.value[j:length(norm.value)])
    phi <- phi.mat[,j]
    phi.prime <- phi.prime.mat[,j]
  } else {
    ## Ignore the initial increasing portion, take the min to the
    ## right of where the initial inflection point occurs
    j <- 1
    while(j < length(norm.stop) && norm.stop[j+1] > norm.stop[j]) j <- j + 1
    j <- j-1 + which.min(norm.stop[j:length(norm.stop)])
    phi <- phi.mat[,j]
    phi.prime <- phi.prime.mat[,j]
  }
  
  console <- printClear(console)
  console <- printPop(console)
  
  if(display.warnings) {
    if(j == iterate.max) warning(" iterate.max reached: increase iterate.max or inspect norm.stop vector")
  }
  
  if(crs.messages) options(crs.messages=FALSE)
  traindata$y <- traindata$y - (fitted(phi.0)-phi)
  
  model <- crs(formula.yz,
               cv="none",
               degree=phi.0$degree,
               segments=phi.0$segments,
               lambda=phi.0$lambda,
               include=phi.0$include,
               kernel=phi.0$kernel,
               basis=phi.0$basis,
               knots=phi.0$knots,
               tau=phi.0$tau,
               deriv=1,
               data=traindata,
               weights=phi.0$weights)
  if(crs.messages) options(crs.messages=TRUE)
  
  model$phi <- phi
  model$phi.prime <- phi.prime
  model$phi.mat <- phi.mat
  model$phi.prime.mat <- phi.prime.mat
  model$num.iterations <- j
  model$norm.stop <- norm.stop
  model$norm.value <- norm.value
  model$convergence <- convergence
  model$starting.values.phi <- starting.values.phi
  model$starting.values.phi.prime <- starting.values.phi.prime
  model$nmulti <- nmulti
  model$ptm <- proc.time() - ptm.start
  
  return(model)
  
}

print.crsivderiv <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

summary.crsivderiv <- function(object, ...) {
  cat("Call:\n")
  print(object$call)
  
  if(!object$kernel) {
    if(is.null(object$tau))
      cat("\nIndicator Bases/B-spline Bases Nonparametric IV Derivative Estimation\n",sep="")
    else
      cat("\nIndicator Bases/B-spline Bases Nonparametric IV Quantile Derivative Estimation\n",sep="")
  } else {
    if(is.null(object$tau))
      cat("\nKernel Weighting/B-spline Bases Nonparametric IV Derivative Estimation\n",sep="")
    else
      cat("\nKernel Weighting/B-spline Bases Nonparametric IV Quantile Derivative Estimation\n",sep="")
  }
  
  if(!is.null(object$tau)) cat(paste("\nQuantile estimated: tau = ",format(object$tau),sep=""),sep="")
  
  cat(paste("\nThere are ",format(object$num.x), " continuous predictors",sep=""),sep="")
  if(!is.null(object$num.z)) cat(paste("\nThere are ",format(object$num.z), " categorical predictors",sep=""),sep="")
  
  for(j in 1:object$num.x)
    cat(paste("\nSpline degree/number of segments for ",format(object$xnames[j]),": ",format(object$degree[j]),"/",format(object$segments[j]),sep=""),sep="")
  if(!is.null(object$include)) for(j in 1:length(object$include))
    cat(paste("\nInclusion indicator for ",format(object$znames[j]),": ",format(object$include[j]),sep=""),sep="")
  if(!is.null(object$lambda)) for(j in 1:length(object$lambda))
    cat(paste("\nBandwidth for ",format(object$znames[j]),": ",format(object$lambda[j]),sep=""),sep="")
  
  cat(paste("\nModel complexity proxy: ", format(object$complexity), sep=""))
  cat(paste("\nKnot type: ", format(object$knots), sep=""))
  if(object$num.x > 1) cat(paste("\nBasis type: ",format(object$basis),sep=""))
  
  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  
  cat(paste("\n\nRegularization method: Landweber-Fridman",sep=""))
  cat(paste("\nNumber of iterations: ", format(object$num.iterations), sep=""))
  cat(paste("\nStopping rule value: ", format(object$norm.stop[length(object$norm.stop)],digits=8), sep=""))
  
  cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))
  cat("\n\n")
}
