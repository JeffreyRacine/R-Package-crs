## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instruments
## xeval: optional evaluation data for the exogenous predictors

## alpha.min: minimum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## alpha.max: maximum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## ... optional arguments for crs()

## This function returns a list with some of the following elements:

## phi: the IV estimator of phi(y)
## alpha:  the Tikhonov regularization parameter
## phi.mat: the matrix with colums phi_1, phi_2 etc. over all iterations
## num.iterations:  the number of Landweber-Fridman iterations
## norm.stop: the vector of values of the objective function used for stopping
## norm.value: the norm not multiplied by the number of iterations
## convergence: a character string indicating whether/why iteration terminated

crsiv <- function(y, ...) UseMethod("crsiv")

crsiv.default <- function(y,
                  z,
                  w,
                  x=NULL,
                  zeval=NULL,
                  weval=NULL,
                  xeval=NULL,
                  alpha=NULL,
                  alpha.max=1.0e-01,
                  alpha.min=1.0e-10,
                  alpha.tol=.Machine$double.eps^0.25,
                  constant=0.5,
                  deriv=0,
                  display.nomad.progress=TRUE,
                  display.warnings=TRUE,
                  iterate.diff.tol=1.0e-08,
                  iterate.max=1000,
                  method=c("Landweber-Fridman","Tikhonov"),
                  opts=list("MAX_BB_EVAL"=10000,
                            "EPSILON"=.Machine$double.eps,
                            "INITIAL_MESH_SIZE"="r1.0e-01",
                            "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                            "MIN_FRAME_SIZE"=paste("r",1,sep=""),
                            "DISPLAY_DEGREE"=0),
                  penalize.iteration=TRUE,
                  smooth.residuals=TRUE,
                  start.from=c("Eyz","EEywz"),
                  starting.values=NULL,
                  stop.on.increase=TRUE,
                  ...) {

  ptm.start <- proc.time()
  old.crs.messages <- getOption("crs.messages")
  on.exit(options(crs.messages = old.crs.messages), add = TRUE)
  crs.messages <- isTRUE(old.crs.messages)
  is.eval.train <- is.null(zeval) && is.null(weval) && is.null(xeval)

  dot.args <- list(...)
  dot.prep <- .crsiv_prepare_dot_args(dot.args)
  weights.arg <- dot.args$weights
  nmulti <- dot.prep$nmulti
  nmulti.loop <- dot.prep$nmulti.loop
  dots.preloop <- dot.prep$dots.preloop
  dots.loop <- dot.prep$dots.loop

  fit.crs <- function(formula, data, dots, degree = NULL, segments = NULL,
                      lambda = NULL, include = NULL, nmulti = NULL) {
    .crsiv_fit_crs(formula = formula,
                   data = data,
                   dots = dots,
                   opts = opts,
                   display.nomad.progress = display.nomad.progress,
                   display.warnings = display.warnings,
                   degree = degree,
                   segments = segments,
                   lambda = lambda,
                   include = include,
                   nmulti = nmulti)
  }

  ## This function was constructed initially by Samuele Centorrino
  ## <samuele.centorrino@univ-tlse1.fr>
  ## the following papers:

  ## A) Econometrica (2011) "Nonparametric Instrumental Regression"
  ## S. Darolles, Y. Fan, J.P. Florens, E. Renault, Volume 79,
  ## 1541-1565.

  ## B) Econometrics Journal (2010), volume 13, pp. S1-S27. doi:
  ## 10.1111/j.1368-423X.2010.00314.x "The practice of non-parametric
  ## estimation by solving inverse problems: the example of
  ## transformation models" Frederique Feve and Jean-Pierre Florens,
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
  ## (Darrolles, Fan, Florens & Renault (2011) it maps W (the
  ## instrument) onto the space of the endogenous regressor Z.

  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.

  ## This function returns TBA (need better error checking!)

  ## phi:   the vector of estimated values for the unknown function at the evaluation points

  tikh <- function(alpha,CZ,CY,Cr.r){
    return(chol2inv(chol(alpha*diag(length(Cr.r)) + CY%*%CZ)) %*% Cr.r) ## This must be computable via ridge... step 1, step 2, same alpha...
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
  ## (Darrolles, Fan, Florens & Renault (2011) it maps W (the
  ## instrument) onto the space of the endogenous regressor Z.

  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.

  ## This function returns TBA (need better error checking!)

  ## phi:   the vector of estimated values for the unknown function at the evaluation points

  ## SSalpha: (scalar) value of the sum of square residuals criterion
  ## which is a function of alpha (see (3.10) of Feve & Florens (2010)

  ## Cr.r is always E.E.y.w.z, r is always E.y.w

  ittik <- function(alpha,CYCZ,Cr.r,r,CZ) {
    invmat <- chol2inv(chol(alpha*diag(length(Cr.r)) + CYCZ))
    tikh.val <- invmat %*% Cr.r
    phi <- tikh.val + alpha * invmat %*% tikh.val
    return(sum((CZ %*% phi - r)^2)/alpha)
  }

  progress.state <- .crs_progress_iv_initialize_state(
    .crs_progress_begin(
      label = .crs_progress_iv_title(),
      domain = "general",
      surface = "iv_solve"
    )
  )
  on.exit({
    progress.state <<- .crs_progress_end(progress.state)
  }, add = TRUE)

  ## Basic error checking

  start.from <- match.arg(start.from)
  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")

  iv_set_stage <- function(label, iteration = NULL) {
    if (!is.null(label)) {
      label <- as.character(label)[1L]
      if (is.na(label) || !nzchar(label)) {
        label <- NULL
      }
    }

    if (!is.null(iteration)) {
      iteration <- suppressWarnings(as.integer(iteration)[1L])
      if (is.na(iteration) || iteration < 1L) {
        iteration <- NULL
      }
    }

    progress.state$iv_object_label <<- label
    progress.state$iv_iteration <<- iteration
    progress.state <<- .crs_progress_step_at(
      state = progress.state,
      now = .crs_progress_now(),
      done = iteration,
      force = TRUE
    )

    invisible(NULL)
  }

  iv_start_label <- function() {
    if (identical(start.from, "Eyz")) "E[y|z]" else "E[E[y|w]|z]"
  }

  iv_residual_stage_label <- function(smooth.residuals) {
    if (smooth.residuals) "E[y-phi(z)|w]" else "E[phi(z)|w]"
  }

  iv_adjoint_stage_label <- function(smooth.residuals) {
    if (smooth.residuals) {
      "E[E[y-phi(z)|w]|z]"
    } else {
      "E[E[y|w]-E[phi(z)|w]|z]"
    }
  }

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(!is.null(x) && NROW(y) != NROW(x)) stop("y and x have differing numbers of rows")
  if(iterate.max < 2) stop("iterate.max must be at least 2")
  if(constant <= 0 || constant >=1) stop("constant must lie in (0,1)")
  if(iterate.diff.tol < 0) stop("iterate.diff.tol must be non-negative")

  ## Cast as data frames

  w <- data.frame(w)
  z <- data.frame(z)
  if(!is.null(x)) x <- data.frame(x)

  ## Check for evaluation data

  if(is.null(zeval)) zeval <- z
  if(is.null(weval)) weval <- w
  if(!is.null(x) && is.null(xeval)) xeval <- x

  method <- match.arg(method)

  if(!is.null(alpha) && alpha <= 0) stop("alpha must be positive")

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

  formula.yw <- as.formula(paste("y ~ ", paste(wnames, collapse= "+")))
  formula.phiw <- as.formula(paste("phi ~ ", paste(wnames, collapse= "+")))
  formula.residw <- as.formula(paste("(y-phi) ~ ", paste(wnames, collapse= "+")))

  if(is.null(x)) {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+")))
    formula.Ephiwz <- as.formula(paste("E.phi.w ~ ", paste(znames, collapse= "+")))
    formula.residwz <- as.formula(paste("residw ~ ", paste(znames, collapse= "+")))
  } else {
    formula.yz <- as.formula(paste("y ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.Eywz <- as.formula(paste("E.y.w ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.Ephiwz <- as.formula(paste("E.phi.w ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
    formula.residwz <- as.formula(paste("residw ~ ", paste(znames, collapse= "+"), " + ", paste(xnames, collapse= "+")))
  }

  if(!is.null(starting.values) && (NROW(starting.values) != NROW(evaldata))) stop(paste("starting.values must be of length",NROW(evaldata)))

  if(method=="Tikhonov") {

    ## Now y=phi(z) + u, hence E(y|w)=E(phi(z)|w) so we need two
    ## bandwidths, one for y on w and one for phi(z) on w (in the
    ## first step we use E(y|w) as a proxy for phi(z) and use
    ## bandwidths for y on w).

    ## convergence flag returned for Landweber-Fridman, not Tikhonov,
    ## but value is required

    convergence <- NULL

    ## First we conduct the regression spline estimator of y on w

    iv_set_stage("E[y|w]")
    .crs_set_messages(crs.messages, FALSE)
    model<-crs(formula.yw,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
    .crs_set_messages(crs.messages, TRUE)

    ## Capture instrument parameters for summary
    degree.w <- model$degree
    segments.w <- model$segments
    lambda.w <- model$lambda
    include.w <- model$include
    num.x.w <- model$num.x
    num.z.w <- model$num.z
    xnames.w <- model$xnames
    znames.w <- model$znames

    E.y.w <- if(is.eval.train) fitted(model) else predict(model,newdata=evaldata,...)
    B <- model.matrix(model$model.lm)
    KYW <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Next, we conduct the regression spline of E(y|w) on z

    iv_set_stage("E[E[y|w]|z]")
    .crs_set_messages(crs.messages, FALSE)
    model <- crs(formula.Eywz,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
    .crs_set_messages(crs.messages, TRUE)
    E.E.y.w.z <- if(is.eval.train) fitted(model) else predict(model,newdata=evaldata,...)
    B <- model.matrix(model$model.lm)
    KYWZ <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Next, we minimize the function ittik to obtain the optimal value
    ## of alpha (here we use the iterated Tikhonov function) to
    ## determine the optimal alpha for the non-iterated scheme. Note
    ## that the function `optimize' accepts bounds on the search (in
    ## this case alpha.min to alpha.max))

    ## E(r|z)=E(E(phi(z)|w)|z)
    ## \phi^\alpha = (\alpha I+CzCw)^{-1}Cr x r

    if(is.null(alpha)) {
      iv_set_stage("alpha")
      alpha <- optimize(ittik, c(alpha.min,alpha.max), tol = alpha.tol, CYCZ = KYWZ %*% KYW, Cr.r = E.E.y.w.z, r = E.y.w, CZ = KYW)$minimum
    }

    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha to get a first stage estimate of phi

    iv_set_stage("phi(z)")
    phi <- as.vector(tikh(alpha, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z))

    ## KYWZ and KZWS no longer used, save memory

    rm(KYW,KYWZ)

    ## Conduct kernel regression of phi(z) on w

    iv_set_stage("E[phi(z)|w]")
    .crs_set_messages(crs.messages, FALSE)
    model <- crs(formula.phiw,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
    .crs_set_messages(crs.messages, TRUE)
    E.phi.w <- if(is.eval.train) fitted(model) else predict(model,newdata=evaldata,...)
    B <- model.matrix(model$model.lm)
    KPHIW <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Conduct kernel regression of E(phi(z)|w) on z

    iv_set_stage("E[E[phi(z)|w]|z]")
    .crs_set_messages(crs.messages, FALSE)
    model <- crs(formula.Ephiwz,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
    .crs_set_messages(crs.messages, TRUE)
    B <- model.matrix(model$model.lm)
    KPHIWZ <- B%*%chol2inv(chol(t(B)%*%B))%*%t(B)

    ## Next, we minimize the function ittik to obtain the optimal value of
    ## alpha (here we use the iterated Tikhonov approach) to determine the
    ## optimal alpha for the non-iterated scheme.

    if(is.null(alpha)) {
      iv_set_stage("alpha")
      alpha <- optimize(ittik,c(alpha.min,alpha.max), tol = alpha.tol, CYCZ = KPHIWZ %*% KPHIW, Cr.r = E.E.y.w.z, r = E.y.w, CZ = KPHIW)$minimum
    }

    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha.

    iv_set_stage("phi(z)")
    phi <- as.vector(tikh(alpha, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.y.w.z))

    iv_set_stage(NULL)

    if(display.warnings) {
      if((alpha-alpha.min)/alpha.min < 0.01) warning(paste(" Tikhonov parameter alpha (",formatC(alpha,digits=4,format="f"),") is close to the search minimum (",alpha.min,")",sep=""))
      if((alpha.max-alpha)/alpha.max < 0.01) warning(paste(" Tikhonov parameter alpha (",formatC(alpha,digits=4,format="f"),") is close to the search maximum (",alpha.max,")",sep=""))
    }

    ## phi.0 is the conditional mean model. We compute lambda =
    ## fitted(phi.0)-phi then transform y via
    ## y.lambda=y-lambda. Here we overwrite y so that we can reuse the
    ## formula. Before that, save the proper residuals and then push
    ## these into the model. June 9 2011 - I am concerned because
    ## phi and the fitted values from this approach are _identical_
    ## (I expected approximately equal).

    ## Feb 21 2012 - JP Florens said the starting point should be
    ## E[E[Y|W]|Z], below we do E[Y|Z]... certainly works, could we
    ## shorten the iterative process?

    .crs_set_messages(crs.messages, FALSE)
    phi.0 <- crs(formula.yz,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)

    residuals.phi <- traindata$y-phi
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
                 deriv=deriv,
                 data=traindata,
                 weights=phi.0$weights)
    .crs_set_messages(crs.messages, TRUE)

    model$residuals <- residuals.phi
    model$phi <- phi
    model$alpha <- alpha
    model$nmulti <- nmulti
    model$ptm <- proc.time() - ptm.start

    ## Attach instrument parameters for summary
    model$degree.w <- degree.w
    model$segments.w <- segments.w
    model$lambda.w <- lambda.w
    model$include.w <- include.w
    model$num.x.w <- num.x.w
    model$num.z.w <- num.z.w
    model$xnames.w <- xnames.w
    model$znames.w <- znames.w

    class(model) <- c("crsiv", "crs")

    return(model)

  } else {

    ## Landweber-Fridman

    ## Create storage vector/matrix

    norm.stop <- numeric()

    ## Compute E(Y|w) for the stopping rule

    iv_set_stage("E[y|w]")

    .crs_set_messages(crs.messages, FALSE)
    model.E.y.w <- crs(formula.yw,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)

    ## Capture instrument parameters for summary
    degree.w <- model.E.y.w$degree
    segments.w <- model.E.y.w$segments
    lambda.w <- model.E.y.w$lambda
    include.w <- model.E.y.w$include
    num.x.w <- model.E.y.w$num.x
    num.z.w <- model.E.y.w$num.z
    xnames.w <- model.E.y.w$xnames
    znames.w <- model.E.y.w$znames

    E.y.w <- if(is.eval.train) fitted(model.E.y.w) else predict(model.E.y.w,newdata=evaldata,...)
    .crs_set_messages(crs.messages, TRUE)

    iv_set_stage(iv_start_label())

    ## Initial value taken from E(E(Y|w)|z) or E(Y|z) or overridden
    ## and passed in, formulae all operate on phi. phi.0.NULL flag set

    .crs_set_messages(crs.messages, FALSE)
    if(is.null(starting.values)) {
      phi.0.NULL <- TRUE
      phi.0 <- crs(formula.yz,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
      ## First compute phi.0 (not passed in) then phi
      if(start.from == "Eyz") {
        ## Start from E(Y|z)
        phi <- if(is.eval.train) fitted(phi.0) else predict(phi.0,newdata=evaldata,...)
      } else {
        ## Start from E(E(Y|w)|z)
        tmp.model <- crs(formula.yw,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
        E.y.w <- fitted(tmp.model)
        model.E.E.y.w.z <- crs(formula.Eywz,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
        phi <- if(is.eval.train) fitted(model.E.E.y.w.z) else predict(model.E.E.y.w.z,newdata=evaldata,...)
      }
    } else {
      phi.0.NULL <- FALSE
      phi.0.input <- starting.values
      ## First compute phi (passed in) then phi.0
      phi <- starting.values
      phi.0 <- crs(formula.yz,opts=opts,data=traindata,display.nomad.progress=display.nomad.progress,display.warnings=display.warnings,...)
    }

    starting.values.phi <- phi

    .crs_set_messages(crs.messages, TRUE)

    iv_set_stage(iv_residual_stage_label(smooth.residuals), iteration = 1L)
    .crs_set_messages(crs.messages, FALSE)
    if(smooth.residuals) {
      traindata$phi <- phi

      model.residw <- fit.crs(formula = formula.residw,
                              data = traindata,
                              dots = dots.preloop)

      ## Capture initial parameters for warm start
      degree.residw <- model.residw$degree
      segments.residw <- model.residw$segments
      lambda.residw <- model.residw$lambda
      include.residw <- model.residw$include

      residw <- if(is.eval.train) fitted(model.residw) else predict(model.residw,newdata=evaldata,...)
      traindata$residw <- residw

      iv_set_stage(iv_adjoint_stage_label(smooth.residuals), iteration = 1L)
      model.predict.residw.z <- fit.crs(formula = formula.residwz,
                                        data = traindata,
                                        dots = dots.preloop)

      ## Capture initial parameters for warm start
      degree.residwz <- model.predict.residw.z$degree
      segments.residwz <- model.predict.residw.z$segments
      lambda.residwz <- model.predict.residw.z$lambda
      include.residwz <- model.predict.residw.z$include

      ## Initialize unused warm start parameters
      degree.phiw <- NULL
      segments.phiw <- NULL
      lambda.phiw <- NULL
      include.phiw <- NULL

    } else {
      traindata$phi <- phi

      model.E.phi.w <- fit.crs(formula = formula.phiw,
                               data = traindata,
                               dots = dots.preloop)

      ## Capture initial parameters for warm start
      degree.phiw <- model.E.phi.w$degree
      segments.phiw <- model.E.phi.w$segments
      lambda.phiw <- model.E.phi.w$lambda
      include.phiw <- model.E.phi.w$include

      residw <- (if(is.eval.train) fitted(model.E.y.w) else predict(model.E.y.w,newdata=evaldata,...)) -
                (if(is.eval.train) fitted(model.E.phi.w) else predict(model.E.phi.w,newdata=evaldata,...))
      traindata$residw <- residw

      iv_set_stage(iv_adjoint_stage_label(smooth.residuals), iteration = 1L)
      model.predict.residw.z <- fit.crs(formula = formula.residwz,
                                        data = traindata,
                                        dots = dots.preloop)

      ## Capture initial parameters for warm start
      degree.residwz <- model.predict.residw.z$degree
      segments.residwz <- model.predict.residw.z$segments
      lambda.residwz <- model.predict.residw.z$lambda
      include.residwz <- model.predict.residw.z$include

      ## Initialize unused warm start parameters
      degree.residw <- NULL
      segments.residw <- NULL
      lambda.residw <- NULL
      include.residw <- NULL

    }
    .crs_set_messages(crs.messages, TRUE)

    if (phi.0.NULL) {
      phi <- (if(is.eval.train) fitted(phi.0) else predict(phi.0,newdata=evaldata,...)) +
             constant*(if(is.eval.train) fitted(model.predict.residw.z) else predict(model.predict.residw.z,newdata=evaldata,...))
    } else {
      phi <- phi.0.input + constant*(if(is.eval.train) fitted(model.predict.residw.z) else predict(model.predict.residw.z,newdata=evaldata,...))
    }

    phi.mat <- matrix(NA, nrow = length(phi), ncol = iterate.max)
    phi.mat[,1] <- phi
    if (!is.null(weights.arg)) {
      weights <- weights.arg
    } else {
      weights <- rep(1, length(y))
    }
    sum_w_Eyw2 <- sum(weights*E.y.w^2)
    norm.stop[1] <- sum(weights*residw^2)/sum_w_Eyw2

    convergence <- "ITERATE_MAX"
    if (iterate.max > 1L) for (j in seq.int(2L, iterate.max)) {

      iv_set_stage(iv_residual_stage_label(smooth.residuals), iteration = j)

      .crs_set_messages(crs.messages, FALSE)
      if(smooth.residuals) {
        traindata$phi <- phi

        model.residw <- fit.crs(formula = formula.residw,
                                data = traindata,
                                dots = dots.loop,
                                degree = degree.residw,
                                segments = segments.residw,
                                lambda = lambda.residw,
                                include = include.residw,
                                nmulti = nmulti.loop)

        degree.residw <- model.residw$degree
        segments.residw <- model.residw$segments
        lambda.residw <- model.residw$lambda
        include.residw <- model.residw$include

        residw <- if(is.eval.train) fitted(model.residw) else predict(model.residw,newdata=evaldata,...)
        traindata$residw <- residw

        iv_set_stage(iv_adjoint_stage_label(smooth.residuals), iteration = j)
        model.predict.residw.z <- fit.crs(formula = formula.residwz,
                                          data = traindata,
                                          dots = dots.loop,
                                          degree = degree.residwz,
                                          segments = segments.residwz,
                                          lambda = lambda.residwz,
                                          include = include.residwz,
                                          nmulti = nmulti.loop)

        degree.residwz <- model.predict.residw.z$degree
        segments.residwz <- model.predict.residw.z$segments
        lambda.residwz <- model.predict.residw.z$lambda
        include.residwz <- model.predict.residw.z$include

      } else {
        traindata$phi <- phi

        model.E.phi.w <- fit.crs(formula = formula.phiw,
                                 data = traindata,
                                 dots = dots.loop,
                                 degree = degree.phiw,
                                 segments = segments.phiw,
                                 lambda = lambda.phiw,
                                 include = include.phiw,
                                 nmulti = nmulti.loop)

        degree.phiw <- model.E.phi.w$degree
        segments.phiw <- model.E.phi.w$segments
        lambda.phiw <- model.E.phi.w$lambda
        include.phiw <- model.E.phi.w$include

        residw <- (if(is.eval.train) fitted(model.E.y.w) else predict(model.E.y.w,newdata=evaldata,...)) -
                  (if(is.eval.train) fitted(model.E.phi.w) else predict(model.E.phi.w,newdata=evaldata,...))
        traindata$residw <- residw

        iv_set_stage(iv_adjoint_stage_label(smooth.residuals), iteration = j)
        model.predict.residw.z <- fit.crs(formula = formula.residwz,
                                          data = traindata,
                                          dots = dots.loop,
                                          degree = degree.residwz,
                                          segments = segments.residwz,
                                          lambda = lambda.residwz,
                                          include = include.residwz,
                                          nmulti = nmulti.loop)

        degree.residwz <- model.predict.residw.z$degree
        segments.residwz <- model.predict.residw.z$segments
        lambda.residwz <- model.predict.residw.z$lambda
        include.residwz <- model.predict.residw.z$include
      }
      .crs_set_messages(crs.messages, TRUE)

      phi <- phi + constant*(if(is.eval.train) fitted(model.predict.residw.z) else predict(model.predict.residw.z,newdata=evaldata,...))
      phi.mat[,j] <- phi

      norm.raw <- sum(weights * residw^2) / sum_w_Eyw2
      norm.stop[j] <- if (penalize.iteration) j * norm.raw else norm.raw

      ## The number of iterations in LF is asymptotically equivalent
      ## to 1/alpha (where alpha is the regularization parameter in
      ## Tikhonov).  Plus the criterion function we use is increasing
      ## for very small number of iterations. So we need a threshold
      ## after which we can pretty much confidently say that the
      ## stopping criterion is decreasing.  In Darolles et al. (2011)
      ## \alpha ~ O(N^(-1/(min(beta,2)+2)), where beta is the so
      ## called qualification of your regularization method. Take the
      ## worst case in which beta = 0 and then the number of
      ## iterations is ~ N^0.5.

      if(j > round(sqrt(nrow(traindata))) && !is.monotone.increasing(norm.stop)) {

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

    }

    phi.mat <- phi.mat[, seq_along(norm.stop), drop = FALSE]
    stop.pick <- .crsiv_select_stop_index(norm.stop)
    norm.value <- stop.pick$norm.value

    ## Extract minimum, and check for monotone increasing function and
    ## issue warning in that case. Otherwise allow for an increasing
    ## then decreasing (and potentially increasing thereafter) portion
    ## of the stopping function, ignore the initial increasing portion,
    ## and take the min from where the initial inflection point occurs
    ## to the length of norm.stop

    if(stop.pick$monotone.failure) {
      .crsiv_warn_monotone_increasing(display.warnings)
      convergence <- "FAILURE_MONOTONE_INCREASING"
    }
    j <- stop.pick$index
    phi <- phi.mat[,j]

    ## phi.0 is the conditional mean model. We compute lambda =
    ## fitted(phi.0)-phi then transform y via
    ## y.lambda=y-lambda. Here we overwrite y so that we can reuse the
    ## formula. Before that, save the proper residuals and then push
    ## these into the model. June 9 2011 - I am concerned because
    ## phi and the fitted values from this approach are _identical_
    ## (I expected approximately equal).

    residuals.phi <- traindata$y-phi
    traindata$y <- traindata$y - (fitted(phi.0)-phi)

    .crs_set_messages(crs.messages, FALSE)
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
                 deriv=deriv,
                 data=traindata,
                 weights=phi.0$weights)
    .crs_set_messages(crs.messages, TRUE)

    model$residuals <- residuals.phi
    model$phi <- phi
    model$phi.mat <- phi.mat
    model$num.iterations <- j
    model$norm.stop <- norm.stop
    model$norm.value <- norm.value
    model$convergence <- convergence
    model$starting.values.phi <- starting.values.phi
    model$nmulti <- nmulti
    model$ptm <- proc.time() - ptm.start

    ## Attach instrument parameters for summary
    model$degree.w <- degree.w
    model$segments.w <- segments.w
    model$lambda.w <- lambda.w
    model$include.w <- include.w
    model$num.x.w <- num.x.w
    model$num.z.w <- num.z.w
    model$xnames.w <- xnames.w
    model$znames.w <- znames.w

    class(model) <- c("crsiv", "crs")

    iv_set_stage(NULL)

    .crsiv_warn_iterate_max(display.warnings, j, iterate.max)

    return(model)

  }

}

print.crsiv <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

summary.crsiv <- function(object, ...) {
  cat("Call:\n")
  print(object$call)

  if(!object$kernel) {
    if(is.null(object$tau))
      cat("\nNonparametric Instrumental Spline Regression\n",sep="")
    else
      cat("\nNonparametric Instrumental Spline Quantile Regression\n",sep="")
  } else {
    if(is.null(object$tau))
      cat("\nNonparametric Instrumental Spline Regression (Kernel Weighting)\n",sep="")
    else
      cat("\nNonparametric Instrumental Spline Quantile Regression (Kernel Weighting)\n",sep="")
  }

  if(!is.null(object$tau)) cat(paste("\nQuantile estimated: tau = ",format(object$tau),sep=""),sep="")

  cat(paste("\nThere are ",format(object$num.x), " continuous predictors",sep=""),sep="")
  if(!is.null(object$num.z)) cat(paste("\nThere are ",format(object$num.z), " categorical predictors",sep=""),sep="")

  for(j in seq_len(object$num.x))
    cat(paste("\nSpline degree/number of segments for ",format(object$xnames[j]),": ",format(object$degree[j]),"/",format(object$segments[j]),sep=""),sep="")
  if(!is.null(object$include)) for(j in seq_along(object$include))
    cat(paste("\nInclusion indicator for ",format(object$znames[j]),": ",format(object$include[j]),sep=""),sep="")
  if(!is.null(object$lambda)) for(j in seq_along(object$lambda))
    cat(paste("\nBandwidth for ",format(object$znames[j]),": ",format(object$lambda[j]),sep=""),sep="")

  if(!is.null(object$num.x.w)) {
    for(j in seq_len(object$num.x.w))
      cat(paste("\nSpline degree/number of segments for ",format(object$xnames.w[j]),": ",format(object$degree.w[j]),"/",format(object$segments.w[j]),sep=""),sep="")
  }
  if(!is.null(object$num.z.w)) {
    if(!is.null(object$include.w)) for(j in seq_along(object$include.w))
      cat(paste("\nInclusion indicator for ",format(object$znames.w[j]),": ",format(object$include.w[j]),sep=""),sep="")
    if(!is.null(object$lambda.w)) for(j in seq_along(object$lambda.w))
      cat(paste("\nBandwidth for ",format(object$znames.w[j]),": ",format(object$lambda.w[j]),sep=""),sep="")
  }

  cat(paste("\nModel complexity proxy: ", format(object$complexity), sep=""))
  cat(paste("\nKnot type: ", format(object$knots), sep=""))
  if(object$num.x > 1) cat(paste("\nBasis type: ",format(object$basis),sep=""))

  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))

  if(!is.null(object$alpha)) {
    cat(paste("\n\nRegularization method: Tikhonov",sep=""))
    cat(paste("\nTikhonov parameter (alpha): ", format(object$alpha,digits=8), sep=""))
  } else {
    cat(paste("\n\nRegularization method: Landweber-Fridman",sep=""))
    cat(paste("\nNumber of iterations: ", format(object$num.iterations), sep=""))
    cat(paste("\nStopping rule value: ", format(object$norm.stop[length(object$norm.stop)],digits=8), sep=""))
  }

  cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))
  cat("\n\n")
}

plot.crsiv <- function(x,
                       plot.data = FALSE,
                       ci = FALSE,
                       deriv = FALSE,
                       xtrim = 0.0,
                       ...) {

  object <- x

  ## We only support univariate endogenous predictor z
  if(object$num.x > 1 || !is.null(object$num.z)) stop(" only univariate z is supported")
  if(!is.logical(plot.data) || length(plot.data) != 1L || is.na(plot.data)) stop(" plot.data must be TRUE/FALSE")
  if(!is.logical(ci) || length(ci) != 1L || is.na(ci)) stop(" ci must be TRUE/FALSE")
  if(!is.logical(deriv) || length(deriv) != 1L || is.na(deriv)) stop(" deriv must be TRUE/FALSE")
  if(!is.numeric(xtrim) || length(xtrim) != 1L || is.na(xtrim)) stop(" xtrim must be a scalar in [0, 0.5)")
  if(xtrim < 0 || xtrim >= 0.5) stop(" xtrim must be in [0, 0.5)")

  z <- object$xz[,1]
  y <- object$y
  phi <- object$phi

  zname <- object$xnames[1]
  yname <- "y" ## Default
  dot_env <- list2env(list(...), parent = emptyenv())

  consume.dot <- function(name, default) {
    if(exists(name, envir = dot_env, inherits = FALSE)) {
      value <- get(name, envir = dot_env, inherits = FALSE)
      rm(list = name, envir = dot_env)
      return(value)
    }
    default
  }

  remaining.dots <- function() {
    as.list.environment(dot_env, all.names = TRUE)
  }

  trim.range <- function(v, trim) {
    if(trim == 0) return(range(v, na.rm = TRUE))
    as.numeric(stats::quantile(v,
                               probs = c(trim, 1 - trim),
                               names = FALSE,
                               na.rm = TRUE,
                               type = 7))
  }

  xlim <- trim.range(z, xtrim)
  keep <- (z >= xlim[1]) & (z <= xlim[2])

  if(deriv) {
    if(is.null(object$deriv.mat)) stop(" deriv.mat not found: was crsiv called with deriv > 0?")
    phi.prime <- object$deriv.mat[keep,1]
    z.plot <- z[keep]
    ord <- order(z.plot)
    z.plot <- z.plot[ord]
    phi.prime <- phi.prime[ord]

    deriv.lwr <- NULL
    deriv.upr <- NULL
    if(ci) {
      if(is.null(object$deriv.mat.lwr) || is.null(object$deriv.mat.upr)) {
        warning("derivative confidence bounds not found; plotting derivative without intervals")
      } else {
        deriv.lwr <- object$deriv.mat.lwr[keep,1][ord]
        deriv.upr <- object$deriv.mat.upr[keep,1][ord]
      }
    }

    ci.lty <- consume.dot("ci.lty", 2)
    ci.lwd <- consume.dot("ci.lwd", 1)
    ci.col <- consume.dot("ci.col", 2)

    plot.args <- c(
      list(x = z.plot,
           y = phi.prime,
           type = consume.dot("type", "l"),
           xlab = consume.dot("xlab", zname),
           ylab = consume.dot("ylab", paste("d", yname, "/d", zname, sep = "")),
           xlim = consume.dot("xlim", xlim),
           ylim = consume.dot("ylim",
                              range(c(phi.prime, deriv.lwr, deriv.upr), na.rm = TRUE)),
           lwd = consume.dot("lwd", 2),
           col = consume.dot("col", 1)),
      remaining.dots()
    )
    do.call(graphics::plot, plot.args)

    if(!is.null(deriv.lwr) && !is.null(deriv.upr)) {
      graphics::lines(z.plot, deriv.lwr,
                      lty = ci.lty,
                      lwd = ci.lwd,
                      col = ci.col)
      graphics::lines(z.plot, deriv.upr,
                      lty = ci.lty,
                      lwd = ci.lwd,
                      col = ci.col)
    }

  } else {
    z.plot <- z[keep]
    phi.plot <- phi[keep]
    ord <- order(z.plot)
    z.plot <- z.plot[ord]
    phi.plot <- phi.plot[ord]

    if(plot.data) {
      y.plot <- y[keep]
      line.lwd <- consume.dot("line.lwd", 2)
      line.lty <- consume.dot("line.lty", 1)
      line.col <- consume.dot("line.col", 1)
      plot.args <- c(
        list(x = z.plot,
             y = y.plot[ord],
             xlab = consume.dot("xlab", zname),
             ylab = consume.dot("ylab", yname),
             xlim = consume.dot("xlim", xlim),
             ylim = consume.dot("ylim", range(c(y.plot, phi.plot), na.rm = TRUE)),
             type = consume.dot("type", "p"),
             col = consume.dot("col", "lightgrey")),
        remaining.dots()
      )
      do.call(graphics::plot, plot.args)
      graphics::lines(z.plot, phi.plot,
                      lwd = line.lwd,
                      lty = line.lty,
                      col = line.col)
    } else {
      plot.args <- c(
        list(x = z.plot,
             y = phi.plot,
             type = consume.dot("type", "l"),
             xlab = consume.dot("xlab", zname),
             ylab = consume.dot("ylab", yname),
             xlim = consume.dot("xlim", xlim),
             ylim = consume.dot("ylim", range(phi.plot, na.rm = TRUE)),
             lwd = consume.dot("lwd", 2),
             col = consume.dot("col", 1)),
        remaining.dots()
      )
      do.call(graphics::plot, plot.args)
    }
  }

}
