## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instruments
## xeval: optional evaluation data for the exogenous predictors

## alpha.min: minimum value when conducting 1-dimensional search for
##            optimal Tihhonov regularization parameter alpha

## alpha.max: maximum value when conducting 1-dimensional search for
##            optimal Tihhonov regularization parameter alpha

## ... optional arguments for crs()

## This function returns a list with the following elements:

## phihat: the IV estimator of phi(y)
## alpha:  the Tikhonov regularization parameter
## num.iterations:  the number of Landweber-Fridman iterations
## norm.stop: the vector of values of the objective function used for stopping

crsiv <- function(y,
                  z,
                  w,
                  x=NULL,
                  zeval=NULL,
                  weval=NULL,
                  xeval=NULL,
                  alpha.min=1.0e-10,
                  alpha.max=1.0e-01,
                  alpha.tol=.Machine$double.eps^0.25,
                  iterate.max=100,
                  iterate.tol=1.0e-05,
                  constant=0.5,
                  method=c("Landweber-Fridman","Tikhonov"),
                  ...) {

  ## This function was constructed initially by Samuele Centorrino
  ## <samuele.centorrino@univ-tlse1.fr>
  ## the following papers:
  
  ## A) Econometrica (forthcoming, article date February 25 2011)
  
  ## "Nonparametric Instrumental Regression"
  ## S. Darolles, Y. Fan, J.P. Florens, E. Renault
  
  ## B) Econometrics Journal (2010), volume 13, pp. S1–S27. doi:
  ## 10.1111/j.1368-423X.2010.00314.x
  
  ## "The practice of non-parametric estimation by solving inverse
  ## problems: the example of transformation models"
  
  ## FREDERIQUE FEVE AND JEAN-PIERRE FLORENS
  ## IDEI and Toulouse School of Economics, Universite de Toulouse
  ## Capitole 21 alle de de Brienne, 31000 Toulouse, France. E-mails:
  ## feve@cict.fr, florens@cict.fr
  
  ## It was modified by Jeffrey S. Racine <racinej@mcmaster.ca> and all
  ## errors remain my responsibility. I am indebted to Samuele and the
  ## Toulouse School of Economics for their generous hospitality.
  
  ## First we require two functions, the first that conducts Regularized
  ## Tikhonov Regression' (aka Ridge Regression)
  
  ## This function conducts regularized Tikhonov regression which
  ## corresponds to (3.9) in Feve & Florens (2010).
  
  ## This function accepts as arguments
  
  ## alpha: penalty
  ## CZ:    row-normalized kernel weights for the `independent' variable
  ## CY:    row-normalized kernel weights for the `dependent' variable
  ## Cr:    row-normalized kernel weights for the `instrument/endogenous' variable (see NOTE below)
  ## r:     vector of conditional expectations (z can be E(Z|z) - see NOTE below)
  
  ## NOTE: for Cr, in the transformation model case treated in Feve &
  ## Florens (2010) this maps Z onto the Y space. In the IV case
  ## (Darrolles, Fan, Florens & Renault (2011, forthcoming Econometrica)
  ## it maps W (the instrument) onto the space of the endogenous
  ## regressor Z.
  
  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.
  
  ## This function returns TBA (need better error checking!)
  
  ## phi:   the vector of estimated values for the unknown function at the evaluation points
  
  tikh <- function(alpha,CZ,CY,Cr.r){
    return(solve(alpha*diag(length(Cr.r)) + CY%*%CZ) %*% Cr.r) ## This must be computable via ridge... step 1, step 2, same alpha...
  }

  ## This function applies the iterated Tikhonov approach which
  ## corresponds to (3.10) in Feve & Florens (2010).
  
  ## This function accepts as arguments
  
  ## alpha: penalty
  ## CZ:    row-normalized kernel weights for the `independent' variable
  ## CY:    row-normalized kernel weights for the `dependent' variable
  ## Cr:    row-normalized kernel weights for the `instrument/endogenous' variable (see NOTE below)
  ## r:     vector of conditional expectations (z can be E(Z|z) - see NOTE below)
  
  ## NOTE: for Cr, in the transformation model case treated in Feve &
  ## Florens (2010) this maps Z onto the Y space. In the IV case
  ## (Darrolles, Fan, Florens & Renault (2011, forthcoming Econometrica)
  ## it maps W (the instrument) onto the space of the endogenous
  ## regressor Z.
  
  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.
  
  ## This function returns TBA (need better error checking!)
  
  ## phi:   the vector of estimated values for the unknown function at the evaluation points
  
  ## SSalpha: (scalar) value of the sum of square residuals criterion
  ## which is a function of alpha (see (3.10) of Feve & Florens (2010)

  ## Cr.r is always E.E.y.w.z, r is always E.y.w
  
  ittik <- function(alpha,CZ,CY,Cr.r,r) {
    invmat <- solve(alpha*diag(length(Cr.r)) + CY%*%CZ)
    tikh.val <- invmat %*% Cr.r
    phi <- tikh.val + alpha * invmat %*% tikh.val ## Not sure about this...
    return(sum((CZ%*%phi - r)^2)/alpha)     ## This is a sum of squared values so CZ%*%phi can be computed with fitted(crs())...
  }

  console <- newLineConsole()

  ## Basic error checking

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(!is.null(x) && NROW(y) != NROW(x)) stop("y and x have differing numbers of rows")
  if(iterate.max < 2) stop("iterate.max must be at least 2")

  ## Check for evaluation data

  if(is.null(zeval)) zeval <- z
  if(is.null(weval)) weval <- w
  if(!is.null(x) && is.null(xeval)) xeval <- x

  method <- match.arg(method)

  ## Set up formulas for multivariate W, Z, and X is provided

  W <- data.frame(w)
  Weval <- data.frame(weval)
  wnames <- paste("w", 1:NCOL(W), sep="")
  names(W) <- wnames
  names(Weval) <- wnames  
  attach(W)

  Z <- data.frame(z)
  Zeval <- data.frame(zeval)  
  znames <- paste("z", 1:NCOL(Z), sep="")
  names(Z) <- znames
  names(Zeval) <- znames  
  attach(Z)

  ## If there exist exogenous regressors X, append these to the
  ## formulas involving Z (can be manually added to W by the user if
  ## desired)

  if(!is.null(x)) {
    X <- data.frame(x)
    Xeval <- data.frame(xeval)    
    xnames <- paste("x", 1:NCOL(X), sep="")
    names(X) <- xnames
    names(Xeval) <- xnames    
    attach(X)
  }

  ## Now create evaluation data

  if(is.null(x)) {
    newdata <- data.frame(Zeval,Weval)
    rm(W,Weval,Z,Zeval)
  } else {
    newdata <- data.frame(Zeval,Weval,Xeval)
    rm(W,Weval,Z,Zeval,X,Xeval)
  }

  formula.yw <- as.formula(paste("y ~ ", paste(wnames, collapse= "+")))
  formula.phihatw <- as.formula(paste("phihat ~ ", paste(wnames, collapse= "+")))  
  formula.residw <- as.formula(paste("(y-phi.j.m.1) ~ ", paste(wnames, collapse= "+")))
  formula.residphi0w <- as.formula(paste("residuals(phi.0) ~ ", paste(wnames, collapse= "+")))    

  if(is.null(x)) {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+")))
    formula.Ephihatwz <- as.formula(paste("E.phihat.w ~ ", paste(znames, collapse= "+")))  
    formula.predictmodelresidphi0z <- as.formula(paste("predict(model.residphi0,newdata=newdata) ~ ", paste(znames, collapse= "+")))
    formula.predictmodelresidwz <- as.formula(paste("predict(model.residw,newdata=newdata) ~ ", paste(znames, collapse= "+")))
  } else {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.Ephihatwz <- as.formula(paste("E.phihat.w ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))  
    formula.predictmodelresidphi0z <- as.formula(paste("predict(model.residphi0,newdata=newdata) ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.predictmodelresidwz <- as.formula(paste("predict(model.residw,newdata=newdata) ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
  }

  if(method=="Tikhonov") {
  
    ## Now y=phi(z) + u, hence E(y|w)=E(phi(z)|w) so we need two
    ## bandwidths, one for y on w and one for phi(z) on w (in the
    ## first step we use E(y|w) as a proxy for phi(z) and use
    ## bandwidths for y on w).
    
    ## First we conduct the regression spline estimator of y on w
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weights and optimal smoothing for E(y|w)...", console)
    model<-crs(formula.yw,...)
    E.y.w <- predict(model,newdata=newdata)
    B <- model.matrix(model$model.lm)
    KYW <- B%*%solve(t(B)%*%B)%*%t(B)
   
    ## Next, we conduct the regression spline of E(y|w) on z
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing weights and optimal smoothing for E(E(y|w)|z)...", console)
    } else {
      console <- printPush("Computing weights and optimal smoothing for E(E(y|w)|z,x)...", console)
    }
    model <- crs(formula.Eywz,...)
    E.E.y.w.z <- predict(model,newdata=newdata)
    B <- model.matrix(model$model.lm)
    KYWZ <- B%*%solve(t(B)%*%B)%*%t(B)
    
    ## Next, we minimize the function ittik to obtain the optimal value
    ## of alpha (here we use the iterated Tikhonov function) to
    ## determine the optimal alpha for the non-iterated scheme. Note
    ## that the function `optimize' accepts bounds on the search (in
    ## this case alpha.min to alpha.max))
    
    ## E(r|z)=E(E(phi(z)|w)|z)
    ## \phi^\alpha = (\alpha I+CzCw)^{-1}Cr x r
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Numerically solving for alpha...", console)
    alpha <- optimize(ittik, c(alpha.min,alpha.max), tol = alpha.tol, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    
    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha to get a first stage estimate of phihat

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing initial phi(z) estimate...", console)
    } else {
      console <- printPush("Computing initial phi(z,x) estimate...", console)
    }
    phihat <- as.vector(tikh(alpha, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z))

    ## KYWZ and KZWS no longer used, save memory
    
    rm(KYW,KYWZ)

    ## Conduct kernel regression of phi(z) on w  
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing optimal smoothing and weights for E(phi(z)|w)...", console)
    } else {
      console <- printPush("Computing optimal smoothing and weights for E(phi(z,x)|w)...", console)
    }
    model <- crs(formula.phihatw,...)
    E.phihat.w <- predict(model,newdata=newdata)
    B <- model.matrix(model$model.lm)
    KPHIW <- B%*%solve(t(B)%*%B)%*%t(B)
    
    ## Conduct kernel regression of E(phi(z)|w) on z
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing optimal smoothing and weights for E(E(phi(z)|w)|z)...", console)
    } else {
      console <- printPush("Computing optimal smoothing and weights for E(E(phi(z,x)|w)|z,x)...", console)
    }
    model <- crs(formula.Ephihatwz,...)
    B <- model.matrix(model$model.lm)
    KPHIWZ <- B%*%solve(t(B)%*%B)%*%t(B)
    
    ## Next, we minimize the function ittik to obtain the optimal value of
    ## alpha (here we use the iterated Tikhonov approach) to determine the
    ## optimal alpha for the non-iterated scheme.
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Iterating and computing the numerical solution for alpha...", console)
    alpha <- optimize(ittik,c(alpha.min,alpha.max), tol = alpha.tol, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    
    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha.

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush("Computing final phi(z) estimate...", console)
    } else {
      console <- printPush("Computing final phi(z,x) estimate...", console)
    }
    phihat <- as.vector(tikh(alpha, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.y.w.z))
    
    console <- printClear(console)
    console <- printPop(console)
    
    if((alpha-alpha.min)/alpha.min < 0.01) warning(paste("Tikhonov parameter alpha (",alpha,") is close to the search minimum (",alpha.min,")",sep=""))
    if((alpha.max-alpha)/alpha.max < 0.01) warning(paste("Tikhonov parameter alpha (",alpha,") is close to the search maximum (",alpha.max,")",sep=""))
    
    return(list(phihat=phihat,alpha=alpha))
    
  } else {
    
    ## Landweber-Fridman

    ## We begin the iteration computing phi.0 and phi.1 directly, then
    ## interate.
    
    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration 1 of maximum ", iterate.max,"...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing and phi(z,x) for iteration 1 of maximum ", iterate.max,"...",sep=""),console)
    }
    phi.0 <- crs(formula.yz,...)
    model.residphi0 <- crs(formula.residphi0w,...)
    model.Eresidphi0.z <- crs(formula.predictmodelresidphi0z,...)
    phi.j.m.1 <- predict(phi.0,newdata=newdata) + predict(model.Eresidphi0.z,newdata=newdata)

    ## For the stopping rule

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing smoothing for stopping rule...",sep=""),console)

    norm.stop <- numeric()
    model.E.y.w <- crs(formula.yw,...)
    E.y.w <- predict(model.E.y.w,newdata=newdata)
    phihat <- phi.j.m.1
    model.E.phi.w <- crs(formula.phihatw,...)
    E.phi.w <- predict(model.E.phi.w,newdata=newdata)
    norm.stop[1] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)

    ascending <- FALSE

    for(j in 2:iterate.max) {

      console <- printClear(console)
      console <- printPop(console)
      if(is.null(x)) {
        console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration ", j, " of maximum ", iterate.max,"...",sep=""),console)
      } else {
        console <- printPush(paste("Computing optimal smoothing and phi(z,x) for iteration ", j, " of maximum ", iterate.max,"...",sep=""),console)
      }

      model.residw <- crs(formula.residw,...)
      model.predict.residw.z <- crs(formula.predictmodelresidwz,...)

      phi.j <- phi.j.m.1 + constant*predict(model.predict.residw.z,newdata=newdata)
      phi.j.m.1 <- phi.j
      phihat <- phi.j

      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing stopping rule for iteration ", j, " of maximum ", iterate.max,"...",sep=""),console)

      ## For the stopping rule (use same smoothing as original)
      model.stop <- crs(formula.phihatw,cv="none",degree=model.E.phi.w$degree,segments=model.E.phi.w$segments,...)
      E.phi.w <- predict(model.stop,newdata=newdata)
      norm.stop[j] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)

      ## If objective increases or we are below stopping tolerance then break

      if((norm.stop[j] > norm.stop[j-1]) || ((norm.stop[j-1]-norm.stop[j]) < iterate.tol)) break()

    }

    console <- printClear(console)
    console <- printPop(console)

    if(j == iterate.max) warning("iterate.max reached: increase iterate.max or inspect norm.stop vector")

    return(list(phihat=phi.j, num.iterations=j, norm.stop=norm.stop))

  }
  
}
