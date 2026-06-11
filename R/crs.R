## This function provides for a complete S3 implementation of
## regression splines with categorical factors using two approaches,
## (i) kernel smoothing, and (ii) indicator function bases. Both
## additive and tensor product bases are supported (default is
## additive, but also see the option basis="auto" that computes both
## and uses that with the smallest cross-validation
## score). Cross-validation (leave-one-out, generalized, and the AIC_c
## method of Hurvich, Simonoff, and Tsai (1998, JRSS B)) can be used
## to select (i) the degree and number of knots (`segments'+1) of the
## basis spline for each continuous predictor, (ii) bandwidth for each
## ordinal/nominal predictor, or (iii) whether or not to include each
## ordinal/nominal predictor's indicator basis. S3 methods include
## fitted, predict, residuals, plot and so forth.

## 2010 (C) Jeffrey S. Racine (racinej@mcmaster.ca).

crs <- function(...) UseMethod("crs")

.crs_default_selector_args <- function() {
  c("cv",
    "cv.df.min",
    "cv.func",
    "cv.threshold",
    "degree.max",
    "degree.min",
    "display.nomad.progress",
    "initial.mesh.size.integer",
    "initial.mesh.size.real",
    "lambda.discrete",
    "lambda.discrete.num",
    "max.bb.eval",
    "min.mesh.size.integer",
    "min.mesh.size.real",
    "min.frame.size.integer",
    "min.frame.size.real",
    "nmulti",
    "opts",
    "random.seed",
    "restarts",
    "segments.max",
    "segments.min",
    "singular.ok")
}

.crs_default_should_select <- function(xz,
                                       kernel,
                                       degree,
                                       segments,
                                       include,
                                       lambda,
                                       dots) {
  xz.df <- as.data.frame(xz)
  dot.names <- names(dots)

  has.numeric <- any(vapply(xz.df, is.numeric, logical(1L)))
  has.categorical <- any(vapply(xz.df, is.factor, logical(1L)))
  has.selector.args <- !is.null(dot.names) &&
    any(nzchar(dot.names) & dot.names %in% .crs_default_selector_args())

  missing.continuous.spec <- has.numeric && (is.null(degree) || is.null(segments))
  missing.categorical.spec <- has.categorical &&
    if (kernel) is.null(lambda) else is.null(include)

  has.selector.args || missing.continuous.spec || missing.categorical.spec
}

.crs_default_delegate <- function(xz,
                                  y,
                                  basis,
                                  complexity,
                                  data.return,
                                  degree,
                                  deriv,
                                  display.nomad.progress,
                                  display.warnings,
                                  include,
                                  kernel,
                                  knots,
                                  lambda,
                                  model.return,
                                  prune,
                                  tau,
                                  weights,
                                  dots) {
  xz.df <- as.data.frame(xz)
  response.name <- make.unique(c(".crs_response", names(xz.df)))[1L]
  data <- xz.df
  data[[response.name]] <- y

  args <- c(
    list(
      formula = reformulate(names(xz.df), response = response.name),
      basis = basis,
      complexity = complexity,
      data = data,
      data.return = data.return,
      degree = degree,
      deriv = deriv,
      display.nomad.progress = display.nomad.progress,
      display.warnings = display.warnings,
      include = include,
      kernel = kernel,
      knots = knots,
      lambda = lambda,
      model.return = model.return,
      prune = prune,
      tau = tau,
      weights = weights
    ),
    dots
  )

  do.call(crs.formula, args)
}

## This function computes the fit and returns the fit, degree
## (vector), segments (vector), and include (vector) for categorical
## predictors. Note that degree of zero and include of zero drop the
## variable from the resulting fit.

crsEst <- function(xz,
                   y,
                   basis=c("auto","additive","tensor","glp"),
                   complexity=c("degree-knots","degree","knots"),
                   data.return=FALSE,
                   degree=NULL,
                   deriv=0,
                   display.nomad.progress=TRUE,
                   display.warnings=TRUE,
                   include=NULL,
                   kernel=TRUE,
                   knots=c("quantiles","uniform","auto"),
                   lambda=NULL,
                   model.return=FALSE,
                   prune=FALSE,
                   prune.index=NULL,
                   segments=NULL,
                   tau=NULL,
                   weights=NULL) {

  ## Take data frame xz and parse into factors (z) and numeric (x).

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)

  if(!kernel) {
    xztmp <- splitFrame(xz)
  } else {
    xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
  }
  x <- xztmp$x
  xnames <- xztmp$xnames
  num.x <- xztmp$num.x
  z <- xztmp$z
  znames <- xztmp$znames
  num.z <- xztmp$num.z
  is.ordered.z <- xztmp$is.ordered.z
  ## The default is kernel==TRUE - this will throw an error with no
  ## categorical predictors so first check
  if(is.null(num.z) && isTRUE(kernel)) kernel <- FALSE
  rm(xztmp)
  if(is.null(z)) {
    include <- NULL
  }

  y <- as.numeric(y)

  ## If weights are provided and there are NA values in the data, we
  ## need only those weights corresponding to the complete cases

  if(!is.null(weights))
    weights <- na.omit(data.frame(xz,y,weights))$weights

  if(!kernel) {

    model <- preditFactorSpline(x=x,
                                y=y,
                                z=z,
                                K=cbind(degree,segments),
                                I=include,
                                knots=knots,
                                basis=basis,
                                prune=prune,
                                tau=tau,
                                weights=weights)

    prune.index <- model$prune.index

    if(deriv > 0) {

      deriv.mat <- xz ## copy for dimension only
      deriv.mat.lwr <- deriv.mat
      deriv.mat.upr <- deriv.mat
      l <- 1 ## num.z
      m <- 1 ## num.x
      for(i in seq_len(ncol(xz))) {
        if(!is.factor(xz[,i])) {
          if(deriv <= degree[m]) {
            tmp <- derivFactorSpline(x=x,
                                     y=y,
                                     z=z,
                                     K=cbind(degree,segments),
                                     I=include,
                                     knots=knots,
                                     basis=basis,
                                     deriv.index=m,
                                     deriv=deriv,
                                     prune.index=prune.index,
                                     tau=tau,
                                     weights=weights)
          } else {
            tmp <- matrix(0,length(y),3)
          }
          deriv.mat[,i] <- tmp[,1]
          deriv.mat.lwr[,i] <- tmp[,2]
          deriv.mat.upr[,i] <- tmp[,3]
          rm(tmp)
          m <- m + 1
        } else {
          ztmp <- z
          ztmp[,l] <- factor(rep(levels(xz[,i])[1],NROW(xz)),levels=levels(xz[,i]),ordered=is.ordered(xz[,i]))

          zpred <- preditFactorSpline(x=x,
                                      y=y,
                                      z=z,
                                      K=cbind(degree,segments),
                                      I=include,
                                      knots=knots,
                                      basis=basis,
                                      prune=prune,
                                      prune.index=prune.index,
                                      tau=tau,
                                      weights=weights)$fitted.values

          zpred.base <- preditFactorSpline(x=x,
                                           y=y,
                                           z=z,
                                           K=cbind(degree,segments),
                                           I=include,
                                           xeval=x,
                                           zeval=ztmp,
                                           knots=knots,
                                           basis=basis,
                                           prune=prune,
                                           prune.index=prune.index,
                                           tau=tau,
                                           weights=weights)$fitted.values

          deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
          deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
          deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

          l <- l + 1
        }
      }

    } else {
      deriv.mat <- NULL
      deriv.mat.lwr <- NULL
      deriv.mat.upr <- NULL
    }

  } else {

    model <- predictKernelSpline(x=x,
                                 y=y,
                                 z=z,
                                 K=cbind(degree,segments),
                                 lambda=lambda,
                                 is.ordered.z=is.ordered.z,
                                 knots=knots,
                                 basis=basis,
                                 model.return=model.return,
                                 tau=tau,
                                 weights=weights)

    prune.index <- NULL

    if(deriv > 0) {

      deriv.mat <- xz ## copy for dimension only
      deriv.mat.lwr <- deriv.mat
      deriv.mat.upr <- deriv.mat
      l <- 1 ## num.z
      m <- 1 ## num.x
      for(i in seq_len(ncol(xz))) {
        if(!is.factor(xz[,i])) {
          if(deriv <= degree[m]) {
            tmp <- derivKernelSpline(x=x,
                                     y=y,
                                     z=z,
                                     K=cbind(degree,segments),
                                     lambda=lambda,
                                     is.ordered.z=is.ordered.z,
                                     knots=knots,
                                     basis=basis,
                                     deriv.index=m,
                                     deriv=deriv,
                                     tau=tau,
                                     weights=weights)
          } else {
            tmp <- matrix(0,length(y),3)
          }
          deriv.mat[,i] <- tmp[,1]
          deriv.mat.lwr[,i] <- tmp[,2]
          deriv.mat.upr[,i] <- tmp[,3]
          rm(tmp)
          m <- m + 1
        } else {

          ztmp <- z
          ztmp[,l] <- rep(sort(unique(z[,l]))[1],NROW(z))

          zpred <- predictKernelSpline(x=x,
                                       y=y,
                                       z=z,
                                       K=cbind(degree,segments),
                                       lambda=lambda,
                                       is.ordered.z=is.ordered.z,
                                       knots=knots,
                                       basis=basis,
                                       model.return=model.return,
                                       tau=tau,
                                       weights=weights)$fitted.values

          zpred.base <- predictKernelSpline(x=x,
                                            y=y,
                                            z=z,
                                            K=cbind(degree,segments),
                                            lambda=lambda,
                                            is.ordered.z=is.ordered.z,
                                            xeval=x,
                                            zeval=ztmp,
                                            knots=knots,
                                            basis=basis,
                                            model.return=model.return,
                                            tau=tau,
                                            weights=weights)$fitted.values

          deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
          deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
          deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

          l <- l + 1
        }
      }

    } else {
      deriv.mat <- NULL
      deriv.mat.lwr <- NULL
      deriv.mat.upr <- NULL
    }

  }

  if(!data.return) {
    x <- NULL
    z <- NULL
  }

  if(is.null(tau)) {
    fitted.values <- model$fitted.values[,1]
    lwr <- model$fitted.values[,2]
    upr <- model$fitted.values[,3]
  } else {
    fitted.values <- model$fitted.values[,1]
    lwr <- NULL
    upr <- NULL
  }

  return(list(fitted.values=fitted.values,
              lwr=lwr,
              upr=upr,
              df.residual=model$df.residual,
              K=cbind(degree,segments),
              degree=degree,
              segments=segments,
              complexity=complexity,
              knots=knots,
              include=include,
              lambda=lambda,
              kernel=kernel,
              basis=basis,
              num.x=num.x,
              num.z=num.z,
              is.ordered.z=is.ordered.z,
              xnames=xnames,
              znames=znames,
              deriv=deriv,
              deriv.mat=deriv.mat,
              deriv.mat.lwr=deriv.mat.lwr,
              deriv.mat.upr=deriv.mat.upr,
              model.lm=model$model,
              hatvalues=model$hatvalues,
              nobs=length(y),
              k=model$rank,
              x=x,
              z=z,
              prune=prune,
              prune.index=prune.index,
              P.hat=model$P.hat,
              tau=tau,
              weights=weights))

}

## Default method - this function takes the minimum arguments (data,
## degree of spline with one element for each column of xz having
## continuous data (presumed default is all xz continuous)).

crs.default <- function(xz,
                        y,
                        basis=c("auto","additive","tensor","glp"),
                        complexity=c("degree-knots","degree","knots"),
                        data.return=FALSE,
                        degree=NULL,
                        deriv=0,
                        display.nomad.progress=TRUE,
                        display.warnings=TRUE,
                        include=NULL,
                        kernel=TRUE,
                        knots=c("quantiles","uniform","auto"),
                        lambda=NULL,
                        model.return=FALSE,
                        prune=FALSE,
                        segments=NULL,
                        tau=NULL,
                        weights=NULL,
                        ...) {
  dots <- list(...)

  if (.crs_default_should_select(xz = xz,
                                 kernel = kernel,
                                 degree = degree,
                                 segments = segments,
                                 include = include,
                                 lambda = lambda,
                                 dots = dots)) {
    est <- .crs_default_delegate(xz = xz,
                                 y = y,
                                 basis = basis,
                                 complexity = complexity,
                                 data.return = data.return,
                                 degree = degree,
                                 deriv = deriv,
                                 display.nomad.progress = display.nomad.progress,
                                 display.warnings = display.warnings,
                                 include = include,
                                 kernel = kernel,
                                 knots = knots,
                                 lambda = lambda,
                                 model.return = model.return,
                                 prune = prune,
                                 tau = tau,
                                 weights = weights,
                                 dots = dots)
    est$call <- match.call(expand.dots = FALSE)
    return(est)
  }

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)

  ## Does the following properly belong here or crsEst?

  est <- crsEst(xz=xz,
                y=y,
                degree=degree,
                segments=segments,
                include=include,
                kernel=kernel,
                lambda=lambda,
                complexity=complexity,
                knots=knots,
                basis=basis,
                deriv=deriv,
                data.return=data.return,
                prune=prune,
                model.return=model.return,
                tau=tau,
                weights=weights,
                display.warnings=display.warnings,
                display.nomad.progress=display.nomad.progress)

  ## Add results to estimated object.

  est$residuals <- y - est$fitted.values
  est$r.squared <- RSQfunc(y,est$fitted.values,weights)
  est$call <- match.call()
  class(est) <- "crs"

  ## Return object of type crs

  return(est)

}

## Here we define the formula and split y (always first column of the
## model frame) from xz (the remaining continuous and
## ordinal/nominal).  nomad/exhaustive search and nmulti controls multiple initial
## points. When nmulti==0, snomadRSolve is used; otherwise
## smultinomadRSolve is used. See ?snomadr
## Jun 4,  2011
##1) degree.max (we have removed  basis.maxdim)
##2) segments.max (we have removed  basis.maxdim)
##3) degree.min (currently 0)
##4) segments.min (currently 1)

crs.formula <- function(formula,
                        basis=c("auto","additive","tensor","glp"),
                        complexity=c("degree-knots","degree","knots"),
                        cv=c("nomad","exhaustive","none"),
                        cv.df.min=1,
                        cv.func=c("cv.ls","cv.gcv","cv.aic"),
                        cv.threshold=1000,
                        data=list(),
                        data.return=FALSE,
                        degree=NULL,
                        degree.max=10,
                        degree.min=0,
                        deriv=0,
                        display.nomad.progress=TRUE,
                        display.warnings=TRUE,
                        include=NULL,
                        initial.mesh.size.integer="1",
                        initial.mesh.size.real="r1.0e-01",
                        kernel=TRUE,
                        knots=c("quantiles","uniform","auto"),
                        lambda=NULL,
                        lambda.discrete=FALSE,
                        lambda.discrete.num=100,
                        max.bb.eval=NULL,
                        min.mesh.size.integer=1,
                        min.mesh.size.real=paste(sqrt(.Machine$double.eps)),
                        min.frame.size.integer=1,
                        min.frame.size.real=1,
                        model.return=FALSE,
                        nmulti=2,
                        opts=list(),
                        prune=FALSE,
                        random.seed=42,
                        restarts=0,
                        segments=NULL,
                        segments.max=10,
                        segments.min=1,
                        singular.ok=FALSE,
                        tau=NULL,
                        weights=NULL,
                        ...) {

  ptm.start <- proc.time()
  opts.supplied <- !missing(opts)
  cv <- match.arg(cv)
  cv.func <- match.arg(cv.func)
  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)

  #  if(!is.null(tau)) {
  #    if(!require(quantreg)) stop(" Error: you must first install the quantreg package")
  #  }

  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  xz <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

  ## Set DISPLAY_DEGREE to 0 if crs.messages=FALSE and DISPLAY_DEGREE
  ## is not provided

  if(!isTRUE(getOption("crs.messages")) && is.null(opts[["DISPLAY_DEGREE"]])) opts$"DISPLAY_DEGREE"=0
  opts.nomad <- if(opts.supplied) opts else list()

  ## If a weights vector is provided and there exists missing data
  ## then the weight vector must be parsed to contain weights
  ## corresponding to the non-missing observations only.

  rows.omit <- as.vector(attr(mf,"na.action"))
  if(!is.null(weights) && !is.null(rows.omit))
    weights <- weights[-rows.omit]

  if(!kernel) {
    xztmp <- splitFrame(xz)
  } else {
    xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
  }
  x <- xztmp$x
  xnames <- xztmp$xnames
  num.x <- xztmp$num.x
  z <- xztmp$z
  znames <- xztmp$znames
  num.z <- xztmp$num.z
  is.ordered.z <- xztmp$is.ordered.z
  ## The default is kernel==TRUE - this will throw an error with no
  ## categorical predictors so first check
  if(is.null(num.z) && isTRUE(kernel)) kernel <- FALSE
  rm(xztmp)
  if(is.null(z)) {
    include <- NULL
  }

  ## Check for dynamic cv and if number of combinations is not overly
  ## large use exhaustive search

  if(cv=="nomad" && is.null(num.z) && ((degree.max-degree.min)*(segments.max-segments.min))**num.x <= cv.threshold) {
    if(display.warnings) warning(" Dynamically changing search from nomad to exhaustive (set cv.threshold=0 to keep NOMAD search)")
    cv <- "exhaustive"
  }

  ## If no degree nor include nor lambda, return cubic spline
  ## (identity bases) or non-smooth model (kernel).

  if(cv.df.min < 1 || cv.df.min > length(y)-1) stop(" cv.df.min must be a positive integer less than n")
  if(!is.null(degree)&&length(degree)!=num.x) stop(" degree vector must be the same length as x")
  if(!is.null(segments)&&length(segments)!=num.x) stop(" segments vector must be the same length as x")
  if(degree.max > 100) stop(paste(" degree.max (",degree.max,") exceeds reasonable value (",100,")",sep=""))
  if(lambda.discrete && (lambda.discrete.num < 1)) stop(" lambda.discrete.num must be a positive integer")

  if(cv=="none"){

    ## When no cross-validation is selected and no defaults are set
    ## for various parameters, we set them to ad hoc defaults and warn
    ## the user to this effect.

    if(is.null(degree) && !is.null(x)) {
      if(display.warnings) warning(paste(" cv=\"none\" selected but no degree provided, using degree=rep(3,num.x): you might consider other degree settings",sep=""),immediate.=TRUE)
      degree <- rep(3,num.x)
    }
    if(is.null(segments) && !is.null(x)) {
      if(display.warnings) warning(paste(" cv=\"none\" selected but no segments provided, using segments=rep(1,num.x): you might consider other segment settings",sep=""),immediate.=TRUE)
      segments <- rep(1,num.x)
    }
    if(is.null(include) && !is.null(z) && !kernel) {
      if(display.warnings) warning(paste(" cv=\"none\" selected but no inclusion for factors indicated, using include=rep(1,num.z): you might consider other include settings",sep=""),immediate.=TRUE)
      include <- rep(1,num.z)
    }
    if(is.null(lambda) && !is.null(z) && kernel) {
      if(display.warnings) warning(paste(" cv=\"none\" selected but no bandwidths for factors indicated, using lambda=rep(0,num.z): you might consider other lambda settings",sep=""),immediate.=TRUE)
      lambda <- rep(0,num.z)
    }

    ## With one continuous predictor all bases are identical, so
    ## simply set the basis to additive and be done (no warning
    ## necessary)

    if(basis == "auto" && num.x == 1) basis <- "additive"

    if(basis == "auto" && num.x > 1) {
      if(display.warnings) warning(paste(" cv=\"none\" selected, basis=\"auto\" changed to basis=\"additive\": you might consider basis=\"tensor\" etc.",sep=""),immediate.=TRUE)
      basis <- "additive"
    }

    if(knots == "auto" && num.x > 1) {
      if(display.warnings) warning(paste(" cv=\"none\" selected, knots=\"auto\" changed to knots=\"quantiles\": you might consider knots=\"uniform\" etc.",sep=""),immediate.=TRUE)
      knots <- "quantiles"
    }

  }

  if(isTRUE(kernel) && isTRUE(prune)) {
    if(display.warnings) warning(" pruning cannot coexist with categorical kernel smoothing (pruning ignored)")
  }
  if(!is.null(tau) && isTRUE(prune)) stop(" pruning is not supported for quantile regression splines")

  ## Check for cv="nomad" and complexity="degree-knots" but
  ## degree.min==degree.max or segments==segments.max

  if((cv=="nomad" && complexity=="degree-knots") && (segments.min==segments.max)) stop("NOMAD search selected with complexity degree-knots but segments.min and segments.max are equal")
  if((cv=="nomad" && complexity=="degree-knots") && (degree.min==degree.max)) stop("NOMAD search selected with complexity degree-knots but degree.min and degree.max are equal")

  ## Check for proper derivative

  if(deriv < 0) stop("derivative order must be a non-negative integer")

  ## Check for logical singular.ok

  if(!is.logical(singular.ok)) stop("singular.ok must be logical (TRUE/FALSE)")

  cv.min <- NULL
  cv.return <- NULL
  cv.maxPenalty <- resolve_cv_maxPenalty(NULL, y, weights = weights, cv.func = cv.func)

  if(!kernel) {

    ## indicator bases and B-spline bases cross-validation

    if(cv=="nomad") {

      max.bb.eval.fr <- if(is.null(max.bb.eval)) 10000 else max.bb.eval

      cv.return <- frscvNOMAD(xz=xz,
                               y=y,
                               degree.max=degree.max,
                               segments.max=segments.max,
                               degree.min=degree.min,
                               segments.min=segments.min,
                               cv.df.min=cv.df.min,
                               complexity=complexity,
                               knots=knots,
                               basis=basis,
                               cv.func=cv.func,
                               degree=degree,
                               segments=segments,
                               random.seed=random.seed,
                               max.bb.eval=max.bb.eval.fr,
                               initial.mesh.size.integer=initial.mesh.size.integer,
                               min.mesh.size.integer=min.mesh.size.integer,
                               min.frame.size.integer=min.frame.size.integer,
                               nmulti=nmulti,
                               opts=opts.nomad,
                               tau=tau,
                               weights=weights,
                               singular.ok=singular.ok,
                               display.nomad.progress=display.nomad.progress,
                               display.warnings=display.warnings)

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min, cv.maxPenalty))) stop(" Search failed: restart with larger nmulti or smaller degree.max  (or degree if provided)")
    }  else if(cv=="exhaustive") {

      cv.return <- frscv(xz=xz,
                          y=y,
                          degree.max=degree.max,
                          segments.max=segments.max,
                          degree.min=degree.min,
                          segments.min=segments.min,
                          complexity=complexity,
                          knots=knots,
                          basis=basis,
                          cv.func=cv.func,
                          degree=degree,
                          segments=segments,
                          tau=tau,
                          weights=weights,
                          singular.ok=singular.ok,
                          display.warnings=display.warnings)

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min, cv.maxPenalty))) stop(" Search failed: restart with smaller degree.max")
    }

  } else {

    ## kernel smooth and B-spline bases cross-validation

    if(cv=="nomad") {

      max.bb.eval.kr <- if(is.null(max.bb.eval)) 1000 else max.bb.eval

      cv.return <- krscvNOMAD(xz=xz,
                               y=y,
                               degree.max=degree.max,
                               segments.max=segments.max,
                               degree.min=degree.min,
                               segments.min=segments.min,
                               cv.df.min=cv.df.min,
                               complexity=complexity,
                               knots=knots,
                               basis=basis,
                               cv.func=cv.func,
                               degree=degree,
                               segments=segments,
                               lambda=lambda,
                               lambda.discrete=lambda.discrete,
                               lambda.discrete.num=lambda.discrete.num,
                               random.seed=random.seed,
                               max.bb.eval=max.bb.eval.kr,
                               initial.mesh.size.real=initial.mesh.size.real,
                               initial.mesh.size.integer=initial.mesh.size.integer,
                               min.mesh.size.real=min.mesh.size.real,
                               min.mesh.size.integer=min.mesh.size.integer,
                               min.frame.size.real=min.frame.size.real,
                               min.frame.size.integer=min.frame.size.integer,
                               nmulti=nmulti,
                               opts=opts.nomad,
                               tau=tau,
                               weights=weights,
                               singular.ok=singular.ok,
                               display.nomad.progress=display.nomad.progress,
                               display.warnings=display.warnings)

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      lambda <- cv.return$lambda
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min, cv.maxPenalty))) stop(" Search failed: restart with larger nmulti or smaller degree.max")

    } else if(cv=="exhaustive") {

      cv.return <- krscv(xz=xz,
                          y=y,
                          degree.max=degree.max,
                          segments.max=segments.max,
                          degree.min=degree.min,
                          segments.min=segments.min,
                          complexity=complexity,
                          knots=knots,
                          basis=basis,
                          cv.func=cv.func,
                          degree=degree,
                          segments=segments,
                          restarts=restarts,
                          tau=tau,
                          weights=weights,
                          singular.ok=singular.ok,
                          display.warnings=display.warnings)

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      lambda <- cv.return$lambda
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min, cv.maxPenalty))) stop(" Search failed: restart with smaller degree.max")

    }

  }

  est <- crs.default(xz=xz,
                      y=y,
                      degree=degree,
                      segments=segments,
                      include=include,
                      kernel=kernel,
                      lambda=lambda,
                      complexity=complexity,
                      knots=knots,
                      basis=basis,
                      deriv=deriv,
                      data.return=data.return,
                      prune=prune,
                      model.return=model.return,
                      tau=tau,
                      weights=weights,
                      display.warnings=display.warnings,
                      ...)


  est$cv.score <- cv.min
  est$call <- match.call()
  est$formula <- formula
  est$terms <- mt
  est$xlevels <- .getXlevels(mt, mf)
  est$xz <- xz
  est$y <- y
  est$prune <- prune
  est$cv.min <- cv.min
  est$cv <- cv
  est$restarts <- restarts
  est$ptm <- proc.time() - ptm.start
  est$nmulti <- nmulti
  if (!is.null(cv.return) && !is.null(cv.return$nomad.restart.contract)) {
    est$nomad.restart.contract <- cv.return$nomad.restart.contract
  }
  if (!is.null(cv.return) && !is.null(cv.return$nomad.best.restart)) {
    est$nomad.best.restart <- cv.return$nomad.best.restart
  }
  if (!is.null(cv.return) && !is.null(cv.return$nomad.restart.objectives)) {
    est$nomad.restart.objectives <- cv.return$nomad.restart.objectives
  }
  if (!is.null(cv.return) && !is.null(cv.return$nomad.restart.evaluations)) {
    est$nomad.restart.evaluations <- cv.return$nomad.restart.evaluations
  }
  est <- .crs_nomad_attach_summary(est, cv.return)

  return(est)

}

## Method for predicting given a new data frame.

predict.crs <- function(object,
                        newdata=NULL,
                        deriv=0,
                        ...) {

  if(is.null(newdata)) {

    ## If no new data provided, return sample fit.
    fitted.values <- fitted(object)
    deriv.mat <- object$deriv.mat

    lwr <- NULL
    upr <- NULL
    deriv.mat.lwr <- NULL
    deriv.mat.upr <- NULL

  } else{

    ## Get training data from object (xz and y) and parse into factors
    ## and numeric.

    basis <- object$basis
    deriv <- object$deriv
    prune <- object$prune
    prune.index <- object$prune.index
    tau <- object$tau
    weights <- object$weights

    xz <- object$xz
    y <- object$y

    ## Divide into factors and numeric

    if(!object$kernel) {
      xztmp <- splitFrame(xz)
    } else {
      xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
    }
    x <- xztmp$x
    z <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z
    rm(xztmp)

    ## Get evaluation data (newdata) and divide into factors and
    ## numeric.

    Terms <- delete.response(terms(object))
    newdata <- model.frame(Terms,newdata,xlev=object$xlevels)

    if(!object$kernel) {
      xztmp <- splitFrame(data.frame(newdata))
    } else {
      xztmp <- splitFrame(data.frame(newdata),factor.to.numeric=TRUE)
    }
    xeval <- xztmp$x
    zeval <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z
    rm(xztmp)

    ## Compute the predicted values.

    if(!object$kernel) {

      ## Get degree vector and include vector.

      complexity <- object$complexity
      knots <- object$knots
      K <- object$K
      degree <- object$degree
      segments <- object$segments
      include <- object$include

      tmp <- preditFactorSpline(x=x,
                                y=y,
                                z=z,
                                K=K,
                                I=include,
                                xeval=xeval,
                                zeval=zeval,
                                basis=basis,
                                knots=knots,
                                prune=prune,
                                prune.index=prune.index,
                                tau=tau,
                                weights=weights)$fitted.values

      fitted.values <- tmp[,1]
      lwr <- tmp[,2]
      upr <- tmp[,3]
      rm(tmp)

      if(deriv > 0) {

        deriv.mat <- matrix(NA,nrow=NROW(newdata),ncol=NCOL(newdata))
        deriv.mat.lwr <- deriv.mat
        deriv.mat.upr <- deriv.mat
        l <- 1 ## num.z
        m <- 1 ## num.x
        for(i in seq_len(ncol(newdata))) {
          if(!is.factor(newdata[,i])) {
            if(deriv <= degree[m]) {
              tmp <- derivFactorSpline(x=x,
                                       y=y,
                                       z=z,
                                       K=K,
                                       I=include,
                                       xeval=xeval,
                                       zeval=zeval,
                                       knots=knots,
                                       basis=basis,
                                       deriv.index=m,
                                       deriv=deriv,
                                       prune.index=prune.index,
                                       tau=tau,
                                       weights=weights)
            } else {
              tmp <- matrix(0,nrow(xeval),3)
            }
            deriv.mat[,i] <- tmp[,1]
            deriv.mat.lwr[,i] <- tmp[,2]
            deriv.mat.upr[,i] <- tmp[,3]
            rm(tmp)
            m <- m + 1
          } else {
            zevaltmp <- zeval
            zevaltmp[,l] <- factor(rep(levels(newdata[,i])[1],NROW(newdata)),levels=levels(newdata[,i]),ordered=is.ordered(newdata[,i]))
            zpred <- preditFactorSpline(x=x,
                                        y=y,
                                        z=z,
                                        K=K,
                                        I=include,
                                        xeval=xeval,
                                        zeval=zeval,
                                        knots=knots,
                                        basis=basis,
                                        prune=prune,
                                        prune.index=prune.index,
                                        tau=tau,
                                        weights=weights)$fitted.values

            zpred.base <- preditFactorSpline(x=x,
                                             y=y,
                                             z=z,
                                             K=K,
                                             I=include,
                                             xeval=xeval,
                                             zeval=zevaltmp,
                                             knots=knots,
                                             basis=basis,
                                             prune=prune,
                                             prune.index=prune.index,
                                             tau=tau,
                                             weights=weights)$fitted.values

            deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
            deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
            deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

            l <- l + 1
          }
        }

      } else {
        deriv.mat <- NULL
        deriv.mat.lwr <- NULL
        deriv.mat.upr <- NULL
      }

    } else {

      ## Get degree vector and lambda vector

      complexity <- object$complexity
      knots <- object$knots
      K <- object$K
      segments <- object$segments
      degree <- object$degree
      lambda <- object$lambda

      is.ordered.z <- object$is.ordered.z

      z <- as.matrix(z)
      zeval <- as.matrix(zeval)

      tmp <- predictKernelSpline(x=x,
                                 y=y,
                                 z=z,
                                 K=K,
                                 lambda=lambda,
                                 is.ordered.z=is.ordered.z,
                                 xeval=xeval,
                                 zeval=zeval,
                                 knots=knots,
                                 basis=basis,
                                 tau=tau,
                                 weights=weights)$fitted.values

      fitted.values <- tmp[,1]
      lwr <- tmp[,2]
      upr <- tmp[,3]
      rm(tmp)

      if(deriv > 0) {

        deriv.mat <- matrix(NA,nrow=NROW(newdata),ncol=NCOL(newdata))
        deriv.mat.lwr <- deriv.mat
        deriv.mat.upr <- deriv.mat
        l <- 1 ## num.z
        m <- 1 ## num.x
        for(i in seq_len(ncol(newdata))) {
          if(!is.factor(newdata[,i])) {
            if(deriv <= degree[m]) {
              tmp <- derivKernelSpline(x=x,
                                       y=y,
                                       z=z,
                                       K=K,
                                       lambda=lambda,
                                       is.ordered.z=is.ordered.z,
                                       xeval=xeval,
                                       zeval=zeval,
                                       knots=knots,
                                       basis=basis,
                                       deriv.index=m,
                                       deriv=deriv,
                                       tau=tau,
                                       weights=weights)
            } else {
              tmp <- matrix(0,nrow(xeval),3)
            }
            deriv.mat[,i] <- tmp[,1]
            deriv.mat.lwr[,i] <- tmp[,2]
            deriv.mat.upr[,i] <- tmp[,3]
            rm(tmp)
            m <- m + 1
          } else {
            zevaltmp <- zeval
            zevaltmp[,l] <- rep(sort(unique(zeval[,l]))[1],NROW(zeval))

            zpred <- predictKernelSpline(x=x,
                                         y=y,
                                         z=z,
                                         K=K,
                                         lambda=lambda,
                                         is.ordered.z=is.ordered.z,
                                         xeval=xeval,
                                         zeval=zeval,
                                         knots=knots,
                                         basis=basis,
                                         tau=tau,
                                         weights=weights)$fitted.values

            zpred.base <- predictKernelSpline(x=x,
                                              y=y,
                                              z=z,
                                              K=K,
                                              lambda=lambda,
                                              is.ordered.z=is.ordered.z,
                                              xeval=xeval,
                                              zeval=zevaltmp,
                                              knots=knots,
                                              basis=basis,
                                              tau=tau,
                                              weights=weights)$fitted.values

            deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
            deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
            deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

            l <- l + 1
          }
        }

      } else {
        deriv.mat <- NULL
        deriv.mat.lwr <- NULL
        deriv.mat.upr <- NULL
      }

    }

  }

  ## Return the predicted values.

  attr(fitted.values, "lwr") <- lwr
  attr(fitted.values, "upr") <- upr
  attr(fitted.values, "deriv.mat") <- deriv.mat
  attr(fitted.values, "deriv.mat.lwr") <- deriv.mat.lwr
  attr(fitted.values, "deriv.mat.upr") <- deriv.mat.upr

  return(fitted.values)

}

## Basic print method.

print.crs <- function(x,
                      ...) {

  cat("Call:\n")
  print(x$call)

}

## print.summary is different from print.

summary.crs <- function(object,
                        sigtest=FALSE,
                        ...) {

  cat("Call:\n")
  print(object$call)
  if(!object$kernel) {
    if(is.null(object$tau))
      cat("\nIndicator Bases/B-spline Bases Regression Spline\n",sep="")
    else
      cat("\nIndicator Bases/B-spline Bases Quantile Regression Spline\n",sep="")
  } else {
    if(is.null(object$tau))
      cat("\nKernel Weighting/B-spline Bases Regression Spline\n",sep="")
    else
      cat("\nKernel Weighting/B-spline Bases Quantile Regression Spline\n",sep="")
  }
  if(!is.null(object$tau)) cat(paste("\nQuantile estimated: tau = ",format(object$tau),sep=""),sep="")
  if(object$num.x==1){
    cat(paste("\nThere is ",format(object$num.x), " continuous predictor",sep=""),sep="")
  } else {
    cat(paste("\nThere are ",format(object$num.x), " continuous predictors",sep=""),sep="")
  }
  if(!is.null(object$num.z)) if(object$num.z==1) {
    cat(paste("\nThere is ",format(object$num.z), " categorical predictor",sep=""),sep="")
  }  else {
    cat(paste("\nThere are ",format(object$num.z), " categorical predictors",sep=""),sep="")
  }
  for(j in seq_len(object$num.x))
    cat(paste("\nSpline degree/number of segments for ",format(object$xnames[j]),": ",format(object$degree[j]),"/",format(object$segments[j]),sep=""),sep="")
  if(!is.null(object$include)) for(j in seq_along(object$include))
    cat(paste("\nInclusion indicator for ",format(object$znames[j]),": ",format(object$include[j]),sep=""),sep="")
  if(!is.null(object$lambda)) for(j in seq_along(object$lambda))
    cat(paste("\nBandwidth for ",format(object$znames[j]),": ",format(object$lambda[j]),sep=""),sep="")
  cat(paste("\nModel complexity proxy: ", format(object$complexity), sep=""))
  cat(paste("\nKnot type: ", format(object$knots), sep=""))
  if(object$num.x > 1) cat(paste("\nBasis type: ",format(object$basis),sep=""))
  if(!object$kernel) cat(paste("\nPruning of final model: ",format(if(object$prune) "TRUE" else "FALSE"),sep=""))
  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  if(is.null(object$tau)) cat(paste("\nRank of model frame: ", format(object$k), sep=""))
  cat(paste("\nTrace of smoother matrix: ", format(round(sum(object$hatvalues))), sep=""))

  if(is.null(object$weights)) {
    if(is.null(object$tau)) cat(paste("\n\nResidual standard error: ", format(sqrt(sum(object$residuals^2)/object$df.residual),digits=4)," on ", format(object$df.residual)," degrees of freedom",sep=""))
    adjusted.r.squared <- 1-(1-object$r.squared)*(length(object$fitted.values)-1)/object$df.residual
    if(is.null(object$tau)) cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4),",   Adjusted R-squared: ",format(adjusted.r.squared,digits=4), sep=""))
    df1 <- round(sum(object$hatvalues))-1
    df2 <- (object$nobs-round(sum(object$hatvalues)))
    F <- (df2/df1)*(sum((object$y-mean(object$y))^2)-sum(residuals(object)^2))/sum(residuals(object)^2)
    if(is.null(object$tau)) cat(paste("\nF-statistic: ", format(F,digits=4), " on ", df1, " and ", df2, " DF, p-value: ", format.pval(pf(F,df1=df1,df2=df2,lower.tail=FALSE),digits=4), sep=""))
    if(!is.null(object$cv.score)) cat(paste("\n\nCross-validation score: ", format(object$cv.score,digits=8), sep=""))
  } else {
    if(is.null(object$tau)) cat(paste("\n\nResidual standard error (weighted): ", format(sqrt(sum((object$residuals^2)*object$weights)/object$df.residual),digits=4)," on ", format(object$df.residual)," degrees of freedom",sep=""))
    adjusted.r.squared <- 1-(1-object$r.squared)*(length(object$fitted.values)-1)/object$df.residual
    if(is.null(object$tau)) cat(paste("\nMultiple R-squared (weighted): ", format(object$r.squared,digits=4),",   Adjusted R-squared (weighted): ",format(adjusted.r.squared,digits=4), sep=""))
    df1 <- round(sum(object$hatvalues))-1
    df2 <- (object$nobs-round(sum(object$hatvalues)))
    F <- (df2/df1)*(sum((object$y-mean(object$y))^2*object$weights)-sum(residuals(object)^2*object$weights))/sum(residuals(object)^2*object$weights)
    if(is.null(object$tau)) cat(paste("\nF-statistic (weighted): ", format(F,digits=4), " on ", df1, " and ", df2, " DF, p-value: ", format.pval(pf(F,df1=df1,df2=df2,lower.tail=FALSE),digits=4), sep=""))
    if(!is.null(object$cv.score)) cat(paste("\n\nCross-validation score (weighted): ", format(object$cv.score,digits=8), sep=""))
  }

  if(object$cv != "none") {
    cat(paste("\nSearch method: ", format(object$cv), sep=""))
    if(identical(object$cv, "nomad"))
      cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
    .crs_nomad_summary_print(object)
  }

  if(sigtest && !object$kernel) {
    cat("\n\nPredictor significance test:\n")
    crs.sigtest(object)
  }

  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))

  cat("\n\n")

}

.crs.SCSrank <- function(x, conf.level = 0.95, alternative = "two.sided") {
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
  DataMatrix <- x
  N <- nrow(DataMatrix)
  k <- round(conf.level * N, 0)
  RankDat <- apply(DataMatrix, 2, rank)
  switch(alternative,
         "two.sided" = {
           W1 <- apply(RankDat, 1, max)
           W2 <- N + 1 - apply(RankDat, 1, min)
           Wmat <- cbind(W1, W2)
           w <- apply(Wmat, 1, max)
           tstar <- round(sort(w)[k], 0)
           SCI <- function(x) {
             sortx <- sort(x)
             cbind(sortx[N + 1 - tstar], sortx[tstar])
           }
           SCS <- t(apply(DataMatrix, 2, SCI))
         },
         "less" = {
           W1 <- apply(RankDat, 1, max)
           tstar <- round(sort(W1)[k], 0)
           SCI <- function(x) {
             sortx <- sort(x)
             cbind(-Inf, sortx[tstar])
           }
           SCS <- t(apply(DataMatrix, 2, SCI))
         },
         "greater" = {
           W2 <- N + 1 - apply(RankDat, 1, min)
           tstar <- round(sort(W2)[k], 0)
           SCI <- function(x) {
             sortx <- sort(x)
             cbind(sortx[N + 1 - tstar], Inf)
           }
           SCS <- t(apply(DataMatrix, 2, SCI))
         })
  colnames(SCS) <- c("lower", "upper")
  list(conf.int = SCS)
}

.crs.bootstrap.bounds <- function(boot.t, alpha, band.type, center) {
  neval <- ncol(boot.t)
  if (band.type == "standard") {
    z <- qnorm(1 - alpha/2)
    se <- sqrt(diag(cov(boot.t)))
    return(cbind(center - z * se, center + z * se))
  }
  if (band.type == "pointwise") {
    return(t(apply(boot.t, 2, quantile, probs = c(alpha/2, 1 - alpha/2))))
  }
  if (band.type == "bonferroni") {
    return(t(apply(boot.t, 2, quantile, probs = c(alpha/(2 * neval), 1 - alpha/(2 * neval)))))
  }
  if (band.type == "simultaneous") {
    return(.crs.SCSrank(boot.t, conf.level = 1 - alpha)$conf.int)
  }
  if (band.type == "all") {
    return(list(
      pointwise = .crs.bootstrap.bounds(boot.t, alpha, "pointwise", center),
      simultaneous = .crs.bootstrap.bounds(boot.t, alpha, "simultaneous", center),
      bonferroni = .crs.bootstrap.bounds(boot.t, alpha, "bonferroni", center)
    ))
  }
  stop("'band.type' must be one of standard, pointwise, bonferroni, simultaneous, all")
}

.crs.draw.all <- function(x, all.bounds, add.legend = TRUE, legend.loc = "topleft") {
  cols <- c(pointwise = "red", simultaneous = "green3", bonferroni = "blue")
  for (bn in c("pointwise", "simultaneous", "bonferroni")) {
    lines(x, all.bounds[[bn]][,1], lty = 2, col = cols[bn])
    lines(x, all.bounds[[bn]][,2], lty = 2, col = cols[bn])
  }
  if (add.legend) {
    legend(legend.loc,
           legend = c("Pointwise", "Simultaneous", "Bonferroni"),
           lty = 2, col = c("red", "green3", "blue"), lwd = 2, bty = "n")
  }
}

.crs.bootstrap.matrix <- function(object,
                                  newdata,
                                  deriv = 0,
                                  deriv.index = 1,
                                  boot.num = 99,
                                  display.warnings = TRUE,
                                  display.nomad.progress = TRUE,
                                  progress.target = NULL) {
  n <- nrow(object$xz)
  prep.activity <- NULL
  if (isTRUE(display.nomad.progress)) {
    prep.activity <- .crs_plot_activity_begin(
      .crs_plot_bootstrap_stage_label(
        stage = "Preparing plot bootstrap",
        target_label = progress.target
      )
    )
    on.exit(.crs_plot_activity_end(prep.activity), add = TRUE)
  }
  pred0 <- predict(object, newdata = newdata, deriv = deriv)
  center <- if (deriv > 0) attr(pred0, "deriv.mat")[,deriv.index] else as.numeric(pred0)
  boot.mat <- matrix(NA_real_, nrow = boot.num, ncol = nrow(newdata))
  progress <- NULL

  if (!is.null(prep.activity)) {
    .crs_plot_activity_end(prep.activity)
    prep.activity <- NULL
  }

  if (isTRUE(display.nomad.progress)) {
    progress <- .crs_plot_stage_progress_begin(
      total = boot.num,
      label = .crs_plot_bootstrap_stage_label(
        stage = "Plot bootstrap",
        target_label = progress.target
      )
    )
    on.exit(.crs_plot_progress_end(progress), add = TRUE)
  }

  for (b in seq_len(boot.num)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    fit.b <- crs.default(
      xz = object$xz[idx,,drop=FALSE],
      y = object$y[idx],
      basis = object$basis,
      complexity = object$complexity,
      degree = object$degree,
      include = object$include,
      kernel = object$kernel,
      knots = object$knots,
      lambda = object$lambda,
      prune = object$prune,
      segments = object$segments,
      tau = object$tau,
      weights = if (is.null(object$weights)) NULL else object$weights[idx],
      display.warnings = display.warnings,
      display.nomad.progress = FALSE
    )
    fit.b$xz <- object$xz[idx,,drop=FALSE]
    fit.b$y <- object$y[idx]
    if (!is.null(object$terms)) fit.b$terms <- object$terms
    if (!is.null(object$xlevels)) fit.b$xlevels <- object$xlevels
    pred.b <- predict(fit.b, newdata = newdata, deriv = deriv)
    boot.mat[b,] <- if (deriv > 0) attr(pred.b, "deriv.mat")[,deriv.index] else as.numeric(pred.b)
    progress <- .crs_plot_progress_tick(progress, done = b, force = (b == 1L))
  }
  list(center = center, boot.mat = boot.mat)
}

.crs_mammen_draws <- function(n, B) {
  a <- (1 - sqrt(5)) / 2
  p.a <- (sqrt(5) + 1) / (2 * sqrt(5))
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1 - a, nrow = n, ncol = B)
  out[u <= p.a] <- a
  out
}

.crs_rademacher_draws <- function(n, B) {
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1, nrow = n, ncol = B)
  out[u <= 0.5] <- -1
  out
}

.crs_plot_normalize_wild <- function(wild = c("rademacher", "mammen")) {
  if(length(wild) > 1L) wild <- wild[1L]
  match.arg(wild, c("mammen", "rademacher"))
}

.crs_wild_chunk_size <- function(n, B) {
  chunk.opt <- getOption("crs.plot.wild.chunk.size")
  if(!is.null(chunk.opt)) {
    chunk.opt <- as.integer(chunk.opt)
    if(length(chunk.opt) != 1L || is.na(chunk.opt) || chunk.opt < 1L) {
      stop("option 'crs.plot.wild.chunk.size' must be a positive integer",
           call. = FALSE)
    }
    return(min(as.integer(B), chunk.opt))
  }
  target.bytes <- 64 * 1024 * 1024
  chunk <- as.integer(floor(target.bytes / (8 * max(1L, as.integer(n)))))
  if(!is.finite(chunk) || is.na(chunk) || chunk < 1L) chunk <- 1L
  min(as.integer(B), chunk)
}

.crs_plot_wild_hat_block_rows <- function(ntrain, neval) {
  target <- getOption("crs.plot.wild.hat.block.bytes", 4 * 1024^2)
  target <- suppressWarnings(as.numeric(target)[1L])
  if(!is.finite(target) || is.na(target) || target <= 0) target <- 4 * 1024^2
  rows <- as.integer(floor(target / (8 * max(1L, as.integer(ntrain)))))
  max(1L, min(as.integer(neval), rows))
}

.crs_plot_wild_dense_hat_enabled <- function(ntrain, neval) {
  threshold <- getOption("crs.plot.wild.dense.hat.threshold.bytes",
                         128 * 1024^2)
  threshold <- suppressWarnings(as.numeric(threshold)[1L])
  if(!is.finite(threshold) || is.na(threshold) || threshold < 0) {
    threshold <- 128 * 1024^2
  }
  as.double(ntrain) * as.double(neval) * 8.0 <= threshold
}

.crs_wild_boot_from_hat <- function(H,
                                    y,
                                    fit.mean,
                                    B,
                                    wild,
                                    display.nomad.progress = TRUE,
                                    progress.label = "Plot bootstrap wild") {
  y <- as.double(y)
  fit.mean <- as.double(fit.mean)
  n <- length(y)
  if(length(fit.mean) != n) {
    stop("length mismatch between fitted means and response for wild bootstrap",
         call. = FALSE)
  }
  B <- as.integer(B)
  if(B < 1L) stop("B must be a positive integer", call. = FALSE)
  wild <- .crs_plot_normalize_wild(wild)
  draw.fun <- if(identical(wild, "mammen")) {
    .crs_mammen_draws
  } else {
    .crs_rademacher_draws
  }
  residuals <- y - fit.mean
  boot.mat <- matrix(NA_real_, nrow = B, ncol = NROW(H))
  chunk.size <- .crs_wild_chunk_size(n = n, B = B)

  progress <- NULL
  if(isTRUE(display.nomad.progress)) {
    progress <- .crs_plot_stage_progress_begin(total = B, label = progress.label)
    on.exit(.crs_plot_progress_end(progress), add = TRUE)
  }

  start <- 1L
  while(start <= B) {
    stopi <- min(B, start + chunk.size - 1L)
    bsz <- stopi - start + 1L
    draws <- draw.fun(n = n, B = bsz)
    ystar <- residuals * draws
    ystar <- ystar + fit.mean
    boot.mat[start:stopi, ] <- t(H %*% ystar)
    progress <- .crs_plot_progress_tick(progress, done = stopi,
                                        force = (start == 1L))
    start <- stopi + 1L
  }

  list(t = boot.mat, t0 = as.vector(H %*% y))
}

.crs.bootstrap.matrix.wild <- function(object,
                                       newdata,
                                       boot.num = 99,
                                       wild = c("rademacher", "mammen"),
                                       display.nomad.progress = TRUE,
                                       progress.target = NULL) {
  if(!is.null(object$tau)) {
    stop("bootstrap=\"wild\" currently supports mean CRS objects only",
         call. = FALSE)
  }
  wild <- .crs_plot_normalize_wild(wild)
  boot.num <- as.integer(boot.num)
  if(boot.num < 1L) stop("B must be a positive integer", call. = FALSE)

  prep.activity <- NULL
  if(isTRUE(display.nomad.progress)) {
    prep.activity <- .crs_plot_activity_begin(
      .crs_plot_bootstrap_stage_label(
        stage = "Preparing plot bootstrap",
        method_label = "wild",
        target_label = progress.target
      )
    )
    on.exit(.crs_plot_activity_end(prep.activity), add = TRUE)
  }

  fit.mean <- as.vector(crshat(object, output = "apply"))
  ntrain <- NROW(object$xz)
  neval <- NROW(newdata)
  center <- as.vector(crshat(object, newdata = newdata, output = "apply"))

  if(!is.null(prep.activity)) {
    .crs_plot_activity_end(prep.activity)
    prep.activity <- NULL
  }

  progress.label <- .crs_plot_bootstrap_stage_label(
    stage = "Plot bootstrap",
    method_label = "wild",
    target_label = progress.target
  )

  if(.crs_plot_wild_dense_hat_enabled(ntrain = ntrain, neval = neval)) {
    H <- crshat(object, newdata = newdata, output = "matrix")
    boot <- .crs_wild_boot_from_hat(
      H = H,
      y = object$y,
      fit.mean = fit.mean,
      B = boot.num,
      wild = wild,
      display.nomad.progress = display.nomad.progress,
      progress.label = progress.label
    )
    return(list(center = center, boot.mat = boot$t))
  }

  boot.mat <- matrix(NA_real_, nrow = boot.num, ncol = neval)
  draw.fun <- if(identical(wild, "mammen")) {
    .crs_mammen_draws
  } else {
    .crs_rademacher_draws
  }
  draws <- draw.fun(n = ntrain, B = boot.num)
  ystar <- (as.double(object$y) - fit.mean) * draws
  ystar <- ystar + fit.mean

  block.rows <- .crs_plot_wild_hat_block_rows(ntrain = ntrain, neval = neval)
  nblocks <- as.integer(ceiling(neval / block.rows))
  progress <- NULL
  if(isTRUE(display.nomad.progress)) {
    progress <- .crs_plot_stage_progress_begin(total = nblocks,
                                               label = progress.label)
    on.exit(.crs_plot_progress_end(progress), add = TRUE)
  }

  start <- 1L
  done <- 0L
  while(start <= neval) {
    stopi <- min(neval, start + block.rows - 1L)
    H <- crshat(object,
                newdata = newdata[start:stopi, , drop = FALSE],
                output = "matrix")
    boot.mat[, start:stopi] <- t(H %*% ystar)
    done <- done + 1L
    progress <- .crs_plot_progress_tick(progress, done = done,
                                        force = (done == 1L))
    start <- stopi + 1L
  }
  list(center = center, boot.mat = boot.mat)
}

.crs_plot_render_surface_rgl <- function(x,
                                         y,
                                         z,
                                         zlim = NULL,
                                         xlab,
                                         ylab,
                                         zlab,
                                         main,
                                         col = NULL,
                                         border = .crs_plot_color("surface_border"),
                                         theta = 45,
                                         phi = 30,
                                         par3d.args = list(),
                                         view3d.args = list(),
                                         persp3d.args = list(),
                                         grid3d.args = list(),
                                         widget.args = list(),
                                         draw.extras = NULL,
                                         data_overlay = FALSE,
                                         data_rug = FALSE,
                                         overlay_x1 = NULL,
                                         overlay_x2 = NULL,
                                         overlay_y = NULL,
                                         display.warnings = TRUE) {
  old.opts <- options(
    rgl.useNULL = TRUE,
    rgl.printRglwidget = TRUE
  )
  on.exit(options(old.opts), add = TRUE)

  old.env <- Sys.getenv("RGL_USE_NULL", unset = NA_character_)
  Sys.setenv(RGL_USE_NULL = "TRUE")
  on.exit({
    if (is.na(old.env)) {
      Sys.unsetenv("RGL_USE_NULL")
    } else {
      Sys.setenv(RGL_USE_NULL = old.env)
    }
  }, add = TRUE)

  if (!isTRUE(suppressWarnings(requireNamespace("rgl", quietly = TRUE)))) {
    if(display.warnings) warning("rgl not installed, option persp.rgl ignored")
    return(invisible(NULL))
  }

  devices.before <- try(rgl::rgl.dev.list(), silent = TRUE)
  if (inherits(devices.before, "try-error") || is.null(devices.before))
    devices.before <- integer(0L)

  opened.dev <- NULL
  cleanup <- function() {
    devices.after <- try(rgl::rgl.dev.list(), silent = TRUE)
    if (!inherits(devices.after, "try-error") && !is.null(devices.after)) {
      new.devices <- setdiff(devices.after, devices.before)
      if (length(new.devices)) {
        for (dev in new.devices) try(rgl::close3d(dev = dev, silent = TRUE), silent = TRUE)
        return(invisible(NULL))
      }
    }
    if (!is.null(opened.dev)) {
      try(rgl::close3d(dev = opened.dev, silent = TRUE), silent = TRUE)
    } else {
      try(rgl::close3d(silent = TRUE), silent = TRUE)
    }
    invisible(NULL)
  }

  tryCatch({
    opened <- rgl::open3d(useNULL = TRUE, silent = TRUE)
    opened.dev <- as.integer(opened[1L])
    on.exit(cleanup(), add = TRUE)

    par3d.call <- .crs_plot_merge_user_args(
      list(windowRect = c(900, 100, 900 + 640, 100 + 640)),
      par3d.args
    )
    do.call(rgl::par3d, par3d.call)

    view3d.call <- .crs_plot_merge_user_args(
      list(theta = theta, phi = phi, fov = 80),
      view3d.args
    )
    do.call(rgl::view3d, view3d.call)

    persp3d.call <- .crs_plot_merge_user_args(
      list(x = x, y = y, z = z,
           zlim = zlim,
           xlab = xlab, ylab = ylab, zlab = zlab,
           ticktype = "detailed",
           border = border,
           color = .crs_plot_rgl_surface_colors(z = z, col = col),
           alpha = 0.6,
           back = "lines",
           main = main),
      persp3d.args
    )
    do.call(rgl::persp3d, persp3d.call)

    grid.side <- c("x", "y+", "z")
    if (!is.null(grid3d.args$side)) {
      grid.side <- grid3d.args$side
      grid3d.args$side <- NULL
    }
    do.call(rgl::grid3d, c(list(grid.side), grid3d.args))

    if (!is.null(draw.extras)) {
      draw.extras()
    }

    if (isTRUE(data_overlay) && !is.null(overlay_x1) &&
        !is.null(overlay_x2) && !is.null(overlay_y)) {
      ok <- is.finite(overlay_x1) & is.finite(overlay_x2) &
        is.finite(overlay_y)
      if (any(ok)) {
        rgl::points3d(overlay_x1[ok], overlay_x2[ok], overlay_y[ok],
                      color = .crs_plot_color("data_overlay", alpha = 1),
                      alpha = 0.35, size = 2)
      }
    }

    if (isTRUE(data_rug) && !is.null(overlay_x1) && !is.null(overlay_x2) &&
        !is.null(zlim)) {
      .crs_plot_draw_floor_rug_rgl(overlay_x1, overlay_x2, zlim)
    }

    widget <- rgl::rglwidget(x = rgl::scene3d())
    if (length(widget.args))
      widget <- do.call(rgl::rglwidget, c(list(x = rgl::scene3d()),
                                         widget.args))
    print(widget)
    invisible(widget)
  }, error = function(e) {
    stop(sprintf("rgl surface renderer failed (%s)", conditionMessage(e)),
         call. = FALSE)
  })
}

plot.crs <- function(x, ...) {
  .crs_plot_regression_1d_public(
    object = x,
    plot.call = match.call(expand.dots = FALSE),
    ...
  )
}

crs.sigtest <- function(object,...) {

  ## This function for the asymptotic significance test can be
  ## airlifted in trivially... trace of the smoother matrix etc. will
  ## deliver correct F stat etc. Left for future 9/1/11 since we have
  ## crssigtest function...

  if(object$kernel) stop(" sigtest is currently available only when kernel=FALSE")

  resolve.call.data <- function(obj) {
    data.expr <- obj$call$data
    if (is.null(data.expr)) {
      return(list())
    }

    env.candidates <- list(
      attr(obj$terms, ".Environment"),
      environment(obj$formula),
      parent.frame()
    )

    for (env in env.candidates) {
      if (!is.environment(env)) {
        next
      }
      data.val <- try(.crs_eval_call(data.expr, env), silent = TRUE)
      if (!inherits(data.val, "try-error")) {
        return(data.val)
      }
    }

    stop("unable to resolve data from fitted call; refit with explicit `data=` in scope")
  }

  model.data <- resolve.call.data(object)
  sg <- list()

  ## Conduct the significance test in order variable by variable

  j.num.x <- 1
  j.num.z <- 1

  for(i in seq_len(NCOL(object$xz))) {

    if(!is.factor(object$xz[,i])) {
      degree <- object$degree
      degree[j.num.x] <- 0
      model.res <- crs(object$formula,cv="none",degree=degree,include=object$include,basis=object$basis,prune=object$prune,data=model.data)
      sg[[i]] <- anova(model.res$model.lm,object$model.lm)
      j.num.x <- j.num.x + 1
    } else {
      include <- object$include
      include[j.num.z] <- 0
      model.res <- crs(object$formula,cv="none",degree=object$degree,include=include,basis=object$basis,prune=object$prune,data=model.data)
      sg[[i]] <- anova(model.res$model.lm,object$model.lm)
      j.num.z <- j.num.z + 1
    }

    cat(paste("Predictor ", format(names(object$xz)[i]), ": Df = ", sg[[i]]$Df[2], ", F = ", format(sg[[i]]$F[2],digits=4), ", Pr(>F) = ", format(sg[[i]][[6]][2],digits=4), "\n", sep=""))

  }

}
