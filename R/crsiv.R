## $Id: crsiv.R,v 1.2 2011/05/26 16:07:04 jracine Exp jracine $

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

## p: order of the local polynomial kernel estimator (p=0 is local
##    constant, p=1 local linear etc.)

## This function returns a list with the following elements:

## phihat: the IV estimator of phi(y)
## alpha:  the Tikhonov regularization parameter

crsiv <- function(y,z,w,yeval=NULL,zeval=NULL,weval=NULL,alpha.min=1.0e-10,alpha.max=1,p=1,tol=.Machine$double.eps^0.25,...) {

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
  
  tikh <- function(alpha,CZ,CY,Cr,r){
    return(solve(alpha*diag(length(r)) + CY%*%CZ) %*% (Cr %*% r))
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
  
  ittik <- function(alpha,CZ,CY,Cr,r) {
    invmat <- solve(alpha*diag(length(r)) + CY%*%CZ)
    phi <- invmat %*% (Cr + alpha*invmat%*%Cr) %*% r
    return((1/alpha)*(crossprod((CZ%*%phi - r),(CZ%*%phi - r))))
  }

  console <- newLineConsole()

  ## Basic error checking

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariat")
  if(NCOL(z) > 1) stop("z must be univariate")
  if(NCOL(w) > 1) stop("w must be univariate")  
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(p < 0) stop("p must be a non-negative integer")

  ## Check for evaluation data

  if(is.null(yeval)) yeval <- y
  if(is.null(zeval)) zeval <- z
  if(is.null(weval)) weval <- w  

  ## First we consider optimal bandwidths for solving E(phi(y)|z)=z
  ## (equation (3.2) in Feve and Florens (2010))
  
  ## Now y=phi(z) + u, hence E(y|w)=E(phi(z)|w) so we need two
  ## bandwidths, one for y on w and one for phi(z) on w (in the first
  ## step we use z on w).
  
  ## NOTE: in what follows E(y|w) is denoted r.
  
  ## First we conduct local polynomial kernel regression of Y on Z to get
  ## the bandwidths.

  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("Computing model and weights for y on w...", console)
  model <- crs(y~w,...)
  B <- model.matrix(model$model.lm)
  KYWs <- B%*%solve(t(B)%*%B)%*%t(B)
  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("Computing model and weights for z on w...", console)
  model <- crs(z~w,...)
  B <- model.matrix(model$model.lm)
  KZWs <- B%*%solve(t(B)%*%B)%*%t(B)
  
  ## r is the conditional expectation of y given w
  
  r <- KYWs%*%y

  ## KYWS no longer used, save memory

  rm(KYWs)
  
  ## define E(r|z)=E(E(phi(z)|w)|z) 
  
  ## We conduct the regression spline of Z on Y we require two
  ## bandwidths, one for r onto z and one for the object in w space
  ## onto z space
  
  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("Computing model and weights for r on z...", console)
  model <- crs(r~z,...)
  B <- model.matrix(model$model.lm)
  KRZs <- B%*%solve(t(B)%*%B)%*%t(B)
  
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
  alpha1 <- optimize(ittik, c(alpha.min,alpha.max), tol = tol, CZ = KZWs, CY = KRZs, Cr = KRZs, r = r)$minimum
  
  ## Finally, we conduct regularized Tikhonov regression using this
  ## optimal alpha.
  
  mized <- as.vector(tikh(alpha1, CZ = KZWs, CY = KRZs, Cr=KRZs, r = r))

  ## KZWS no longer used, save memory

  rm(KZWs)
  
  ## Here we need to generate the kernel weights for the local
  ## polynomial estimator. This corresponds to C_z below (3.8) in Feve &
  ## Florens (2010).
  
  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("Computing model and weights for phi on w...", console)
  model <- crs(mized~w,...)
  B <- model.matrix(model$model.lm)
  KPHWs <- B%*%solve(t(B)%*%B)%*%t(B)
  
  ## Conduct kernel regression of E(phi(z)|w) on z
  
  ## Here we need to generate the kernel weights for the local
  ## polynomial estimator. This corresponds to C_y below (3.8) in Feve &
  ## Florens (2010).
  
  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("Iterating and recomputing model and weights for phi on z...", console)
  model <- crs(as.vector(KPHWs%*%mized)~z,...)
  B <- model.matrix(model$model.lm)
  KPHZs <- B%*%solve(t(B)%*%B)%*%t(B)
  
  ## Next, we minimize the function ittik to obtain the optimal value of
  ## alpha (here we use the iterated Tikhonov approach) to determine the
  ## optimal alpha for the non-iterated scheme.
  
  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("Iterating and recomputing the numerical solution for alpha...", console)
  alpha2 <- optimize(ittik,c(alpha.min,alpha.max), tol = tol, CZ = KPHWs, CY = KPHZs, Cr = KRZs, r = r)$minimum
  
  ## Finally, we conduct regularized Tikhonov regression using this
  ## optimal alpha.

  mized2 <- as.vector(tikh(alpha2, CZ = KPHWs, CY = KPHZs, Cr = KRZs, r = r))
  
  console <- printClear(console)
  console <- printPop(console)

  return(list(phihat=mized2,alpha=alpha2))
  
}


