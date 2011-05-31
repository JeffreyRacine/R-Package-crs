## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictor
## w: instrument

## yeval: optional evaluation data for the univariate outcome
## zeval: optional evaluation data for the endogenous predictor
## weval: optional evaluation data for the instrument

## alpha.min: minimum value when conducting 1-dimensional search for
##            optimal Tihhonov regularization parameter alpha

## alpha.max: maximum value when conducting 1-dimensional search for
##            optimal Tihhonov regularization parameter alpha

## ... optional arguments for crs()

## This function returns a list with the following elements:

## phihat: the IV estimator of phi(y)
## alpha:  the Tikhonov regularization parameter

crsiv <- function(y,
                  z,
                  w,
                  yeval=NULL,
                  zeval=NULL,
                  weval=NULL,
                  alpha.min=1.0e-10,
                  alpha.max=1,
                  tol=.Machine$double.eps^0.25,
                  num.iterations=5,
                  max.iterations=25,
                  constant=0.5,
                  method=c("Landweber-Fridman","Tikhonov"),
                  ...) {

  ## This function was constructed initially by Samuele Centorrino
  ## <samuele.centorrino@univ-tlse1.fr> to reproduce illustrations in
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
    return(solve(alpha*diag(length(Cr.r)) + CY%*%CZ) %*% Cr.r)
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
  
  ittik <- function(alpha,CZ,CY,Cr.r,r) {
    invmat <- solve(alpha*diag(length(Cr.r)) + CY%*%CZ)
    phi <- invmat %*% Cr.r + alpha * invmat %*% invmat %*% Cr.r        
    return(sum((CZ%*%phi - r)^2)/alpha)    
  }

  console <- newLineConsole()

  ## Basic error checking

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariat")
  if(NCOL(z) > 1) stop("z must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(num.iterations < 2) stop("num.iterations must be at least 2")

  ## Check for evaluation data

  if(is.null(yeval)) yeval <- y
  if(is.null(zeval)) zeval <- z
  if(is.null(weval)) weval <- w

  method <- match.arg(method)

  ## Set up formulas for multivariate W

  W <- data.frame(w)
  wnames <- paste("w", 1:NCOL(W), sep="")
  names(W) <- wnames
  attach(W)
  rm(W)
  formula.zw <- as.formula(paste("z ~ ", paste(wnames, collapse= "+")))
  formula.yw <- as.formula(paste("y ~ ", paste(wnames, collapse= "+")))
  formula.phihatw <- as.formula(paste("phihat ~ ", paste(wnames, collapse= "+")))  
  formula.residw <- as.formula(paste("(y-phi.j.m.1) ~ ", paste(wnames, collapse= "+")))
  formula.residphi.0 <- as.formula(paste("residuals(phi.0) ~ ", paste(wnames, collapse= "+")))    

  if(method=="Tikhonov") {
  
    ## Now y=phi(z) + u, hence E(y|w)=E(phi(z)|w) so we need two
    ## bandwidths, one for y on w and one for phi(z) on w (in the first
    ## step we use z on w).
    
    ## First we conduct the regression spline estimator of y on w
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing E(y|w)...", console)
    E.y.w <- fitted(crs(formula.yw,...))
    
    ## Next, we conduct the regression spline of E(y|w) on z
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing model and weights for E(E(y|w)|z)...", console)
    model <- crs(E.y.w~z,...)
    E.E.y.w.z <- fitted(model)
    B <- model.matrix(model$model.lm)
    KRZs <- B%*%solve(t(B)%*%B)%*%t(B)
    
    ## Next, weights for E(z|w)
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing model and weights for E(z|w) (first stage treat z as phi(z))...", console)
    model <- crs(formula.zw,...)
    B <- model.matrix(model$model.lm)
    KZWs <- B%*%solve(t(B)%*%B)%*%t(B)
    
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
    alpha1 <- optimize(ittik, c(alpha.min,alpha.max), tol = tol, CZ = KZWs, CY = KRZs, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    
    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha to get a first stage estimate of phihat
    
    phihat <- as.vector(tikh(alpha1, CZ = KZWs, CY = KRZs, Cr.r = E.E.y.w.z))
    
    ## KZWS no longer used, save memory
    
    rm(KZWs)
    
    ## Conduct kernel regression of phi(z) on w  
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing model and weights for E(phi(z)|w)...", console)
    model <- crs(formula.phihatw,...)
    E.phiyat.w <- fitted(model)
    B <- model.matrix(model$model.lm)
    KPHIWs <- B%*%solve(t(B)%*%B)%*%t(B)
    
    ## Conduct kernel regression of E(phi(z)|w) on z
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Iterating and recomputing model and weights for E(E(phi(z)|w)|z)...", console)
    model <- crs(E.phiyat.w~z,...)
    B <- model.matrix(model$model.lm)
    KPHIZs <- B%*%solve(t(B)%*%B)%*%t(B)
    
    ## Next, we minimize the function ittik to obtain the optimal value of
    ## alpha (here we use the iterated Tikhonov approach) to determine the
    ## optimal alpha for the non-iterated scheme.
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Iterating and recomputing the numerical solution for alpha...", console)
    alpha2 <- optimize(ittik,c(alpha.min,alpha.max), tol = tol, CZ = KPHIWs, CY = KPHIZs, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    
    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha.
    
    phihat2 <- as.vector(tikh(alpha2, CZ = KPHIWs, CY = KPHIZs, Cr.r = E.E.y.w.z))
    
    console <- printClear(console)
    console <- printPop(console)
    
    return(list(phihat=phihat2,alpha=alpha2))
    
  } else {
    
    ## Landweber-Fridman

    ## We begin the iteration computing phi.0 and phi.1 directly, then
    ## interate.
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration 1 of at least ", num.iterations,"...",sep=""),console)
    phi.0 <- crs(y~z,...)
    phi.j.m.1 <- fitted(phi.0) + fitted(crs(fitted(crs(formula.residphi.0,...))~z,...))

    ## For the stopping rule

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing smoothing for stopping rule...",sep=""),console)

    norm.stop <- numeric()
    model.E.y.w <- crs(formula.yw,...)
    E.y.w <- fitted(model.E.y.w)
    phihat <- phi.j.m.1
    model.E.phi.w <- crs(formula.phihatw,...)
    E.phi.w <- fitted(model.E.phi.w)
    norm.stop[1] <- sum((E.y.w-E.phi.w)^2)

    for(j in 2:num.iterations) {
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration ", j, " of at least ", num.iterations,"...",sep=""),console)

      phi.j <- phi.j.m.1 + constant*fitted(crs(fitted(crs(formula.residw,...))~z,...))
      phi.j.m.1 <- phi.j
      phihat <- phi.j

      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing stopping rule for iteration ", j, " of at least ", num.iterations,"...",sep=""),console)

      ## For the stopping rule (use same smoothing as original)
      E.phi.w <- fitted(crs(formula.phihatw,cv="none",degree=model.E.phi.w$degree,segments=model.E.phi.w$segments,...))
      norm.stop[j] <- sum((E.y.w-E.phi.w)^2)
    }

    ## If the last num.iterations normed differences are unchanged, stop

    while((sum(norm.stop[(j-num.iterations+2):j]-norm.stop[(j-num.iterations+1):(j-1)]) != 0) && j < max.iterations) {
      j <- j+1
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing optimal smoothing and phi(z) for iteration ", j, "of a maximum of ", max.iterations, "...",sep=""),console)

      phi.j <- phi.j.m.1 + constant*fitted(crs(fitted(crs(formula.residw,...))~z,...))
      phi.j.m.1 <- phi.j
      phihat <- phi.j

      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing stopping rule for iteration ", j, "of a maximum of ", max.iterations, "...",sep=""),console)

      ## For the stopping rule (use same smoothing as original)
      E.phi.w <- fitted(crs(formula.phihatw,cv="none",degree=model.E.phi.w$degree,segments=model.E.phi.w$segments,...))
      norm.stop[j] <- sum((E.y.w-E.phi.w)^2)
    }
    
    console <- printClear(console)
    console <- printPop(console)

    if(j == max.iterations) warning("max.iterations reached: increase max.iterations or inspect norm.stop vector")

    return(list(phihat=phi.j, num.iterations=j, norm.stop=norm.stop))

  }
  
}
