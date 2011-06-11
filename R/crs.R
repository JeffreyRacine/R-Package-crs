## 3/6/2010, (C) Jeffrey S. Racine (racinej@mcmaster.ca).

## This function provides for a complete S3 implementation of
##  regression splines with categorical factors using two approaches,
##  (i) kernel smoothing, and (ii) indicator function
##  bases. Cross-validation (leave-one-out) can be used to select (i)
##  the degree of the basis spline for each continuous predictor, (ii)
##  bandwidth for each ordinal/nominal predictor, or (iii) whether or
##  not to include each ordinal/nominal predictor's indicator basis.

crs <- function(...) UseMethod("crs")

## This function computes the fit and returns the fit, degree
## (vector), and include (vector) for categorical predictors. Note
## that degree of zero and include of zero drop the variable from the
## resulting fit.

crsEst <- function(xz,
                   y,
                   degree=NULL,
                   segments=NULL,
                   include=NULL,
                   kernel=TRUE,
                   lambda=NULL,
                   kernel.type=c("nominal","ordinal"),
                   complexity=c("degree-knots","degree","knots"),
                   knots=c("quantiles","uniform"),
                   basis=c("additive","tensor","auto"),
                   deriv=0,
                   data.return=FALSE,
                   prune=FALSE,
                   prune.index=NULL) {
  
  ## Take data frame xz and parse into factors (z) and numeric (x).

  kernel.type <- match.arg(kernel.type)
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
  ## The default is kernel==TRUE - this will throw an error with no
  ## categorical predictors so first check
  if(is.null(num.z) && kernel==TRUE) kernel <- FALSE
  rm(xztmp)
  if(is.null(z)) {
    include <- NULL
  }

  y <- as.numeric(y)

  if(!kernel) {

    model <- predict.factor.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=cbind(degree,segments),
                                   I=include,
                                   knots=knots,
                                   basis=basis,
                                   prune=prune)

    prune.index <- model$prune.index

    if(deriv!=0) {

      deriv.mat <- xz ## copy for dimension only
      deriv.mat.lwr <- deriv.mat
      deriv.mat.upr <- deriv.mat
      l <- 1 ## num.z
      m <- 1 ## num.x
      for(i in 1:ncol(xz)) {
        if(!is.factor(xz[,i])) {
          tmp <- deriv.factor.spline(x=x,
                                     y=y,
                                     z=z,
                                     K=cbind(degree,segments),
                                     I=include,
                                     knots=knots,
                                     basis=basis,
                                     deriv.index=m,
                                     deriv=deriv,
                                     prune.index=prune.index)
          deriv.mat[,i] <- tmp[,1]
          deriv.mat.lwr[,i] <- tmp[,2]
          deriv.mat.upr[,i] <- tmp[,3]
          rm(tmp)
          m <- m + 1
        } else {
          ztmp <- z
          ztmp[,l] <- factor(rep(levels(xz[,i])[1],NROW(xz)),levels=levels(xz[,i]))

          zpred <- predict.factor.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=cbind(degree,segments),
                                         I=include,
                                         knots=knots,
                                         basis=basis,
                                         prune=prune,
                                         prune.index=prune.index)$fitted.values

          zpred.base <- predict.factor.spline(x=x,
                                              y=y,
                                              z=z,
                                              K=cbind(degree,segments),
                                              I=include,
                                              xeval=x,
                                              zeval=ztmp,
                                              knots=knots,
                                              basis=basis,
                                              prune=prune,
                                              prune.index=prune.index)$fitted.values

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

    model <- predict.kernel.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=cbind(degree,segments),
                                   lambda=lambda,
                                   kernel.type=kernel.type,
                                   knots=knots,
                                   basis=basis)

    prune.index <- NULL

    if(deriv!=0) {

      deriv.mat <- xz ## copy for dimension only
      deriv.mat.lwr <- deriv.mat
      deriv.mat.upr <- deriv.mat
      l <- 1 ## num.z
      m <- 1 ## num.x
      for(i in 1:ncol(xz)) {
        if(!is.factor(xz[,i])) {
          tmp <- deriv.kernel.spline(x=x,
                                     y=y,
                                     z=z,
                                     K=cbind(degree,segments),
                                     lambda=lambda,
                                     kernel.type=kernel.type,
                                     knots=knots,
                                     basis=basis,
                                     deriv.index=m,
                                     deriv=deriv)

          deriv.mat[,i] <- tmp[,1]
          deriv.mat.lwr[,i] <- tmp[,2]
          deriv.mat.upr[,i] <- tmp[,3]
          rm(tmp)
          m <- m + 1
        } else {
          ztmp <- z
          ztmp[,l] <- as.numeric(factor(rep(levels(xz[,i])[1],NROW(xz)),levels=levels(xz[,i])))

          zpred <- predict.kernel.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=cbind(degree,segments),
                                         lambda=lambda,
                                         kernel.type=kernel.type,
                                         knots=knots,
                                         basis=basis)$fitted.values

          zpred.base <- predict.kernel.spline(x=x,
                                              y=y,
                                              z=z,
                                              K=cbind(degree,segments),
                                              lambda=lambda,
                                              kernel.type=kernel.type,
                                              xeval=x,
                                              zeval=ztmp,
                                              knots=knots,
                                              basis=basis)$fitted.values

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

  return(list(fitted.values=model$fitted.values[,1],
              lwr=model$fitted.values[,2],
              upr=model$fitted.values[,3],
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
              cv.score=mean((y-model$fitted.values[,1])^2/(1-model$hatvalues)^2),
              x=x,
              z=z,
              prune=prune,
              prune.index=prune.index))

}

## Default method - this function takes the minimum arguments (data,
## degree of spline with one element for each column of xz having
## continuous data (presumed default is all xz continuous)).

crs.default <- function(xz,
                        y,
                        degree=NULL,
                        segments=NULL,
                        include=NULL,
                        kernel=TRUE,
                        lambda=NULL,
                        kernel.type=c("nominal","ordinal"),
                        complexity=c("degree-knots","degree","knots"),
                        knots=c("quantiles","uniform"),
                        basis=c("additive","tensor","auto"),
                        deriv=0,
                        data.return=FALSE,
                        prune=FALSE,
                        ...) {

  kernel.type <- match.arg(kernel.type)
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
                kernel.type=kernel.type,
                complexity=complexity,
                knots=knots,
                basis=basis,
                deriv=deriv,
                data.return=data.return,
                prune=prune)

  ## Add results to estimated object.

  est$residuals <- y - est$fitted.values
  est$r.squared <- RSQfunc(y,est$fitted.values)
  est$call <- match.call()
  class(est) <- "crs"

  ## Return object of type crs

  return(est)

}

## Here we define the formula and split y (always first column of the
## model frame) from xz (the remaining continuous and
## ordinal/nominal).  nomad::FALSE exhaustive search nmulti is the
## number for multiple initial points.  if it is bigger than 1, when
## nomad is true, it will call snomadRSolve, otherwise, it will call
## smultinomadRSolve See ?snomadr
# Jun 4,  2011
#1) degree.max (we have removed  basis.maxdim)
#2) segments.max (we have removed  basis.maxdim)
#3) degree.min (currently 0)
#4) segments.min (currently 1)

crs.formula <- function(formula,
                        data=list(),
                        degree=NULL,
                        segments=NULL,
                        include=NULL,
												degree.max=10, 
												segments.max=10, 
												degree.min=0, 
												segments.min=1, 
                        cv=c("nomad","exhaustive","none"),
                        cv.func=c("cv.ls","cv.gcv","cv.aic"),
                        kernel=TRUE,
                        lambda=NULL,
                        kernel.type=c("nominal","ordinal"),
                        complexity=c("degree-knots","degree","knots"),
                        knots=c("quantiles","uniform"),
                        basis=c("additive","tensor","auto"),
                        deriv=0,
                        data.return=FALSE,
                        prune=FALSE,
												restarts=0,
												opts=list("MAX_BB_EVAL"=500,"MIN_MESH_SIZE"="r1.0e-10","INITIAL_MESH_SIZE"="r1.0e-00","MIN_POLL_SIZE"="r1.0e-10"),
												nmulti=0, 
												...) {

  cv <- match.arg(cv)  
  cv.func <- match.arg(cv.func)
  kernel.type <- match.arg(kernel.type)
  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)

  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  y <- model.response(mf)
#  xz <- data.frame(mf[,-1,drop=FALSE])
  ## May 20 2011, trying alternative way for factors...
  xz <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]
#  names(xz) <- names(mf)[-1] ## Case of one predictor has names clobbered

  ### cv needs this? 9/12/2010

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
  ## The default is kernel==TRUE - this will throw an error with no
  ## categorical predictors so first check
  if(is.null(num.z) && kernel==TRUE) kernel <- FALSE
  rm(xztmp)
  if(is.null(z)) {
    include <- NULL
  }

  ## If no degree nor include nor lambda, return cubic spline
  ## (identity bases) or non-smooth model (kernel).

  if(!is.null(degree)&&length(degree)!=num.x) stop(" degree vector must be the same length as x")
  if(!is.null(segments)&&length(segments)!=num.x) stop(" segments vector must be the same length as x")    

	if(cv=="none"){
			if(is.null(degree)&!is.null(x)) degree <- rep(3,num.x)
			if(is.null(segments)&!is.null(x)) segments <- rep(1,num.x)
			if(is.null(include)&!is.null(z)&!kernel) include <- rep(1,num.z)
			if(is.null(lambda)&!is.null(z)&kernel) lambda <- rep(0,num.z)
	}


  if(cv!="none"&&basis!="auto"&&NCOL(xz)>1) warning(paste(" cv specified but basis is ", basis, ": you might consider basis=\"auto\"",sep=""))

  if(kernel==TRUE&&prune==TRUE) warning(" pruning cannot coexist with categorical kernel smoothing (pruning ignored)")

  cv.min <- NULL

  if(!kernel) {

    ## indicator bases and B-spline bases cross-validation
    
    if(cv=="nomad") {
      
      cv.return <- frscvNOMAD(xz=xz,
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
															include=include, 
															opts=opts,
															nmulti=nmulti)

	cv.min <- cv.return$cv.min
	degree <- cv.return$degree
	segments <- cv.return$segments
	include <- cv.return$I
	basis <- cv.return$basis

		}	else if(cv=="exhaustive") {

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
												 segments=segments)


	cv.min <- cv.return$cv.min
	degree <- cv.return$degree
	segments <- cv.return$segments
	include <- cv.return$I
	basis <- cv.return$basis

		}

  } else {

    ## kernel smooth and B-spline bases cross-validation
    
    if(cv=="nomad") {

      cv.return <- krscvNOMAD(xz=xz,
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
															lambda=lambda, 
                              restarts=restarts, 
                              opts=opts,
                              nmulti=nmulti)
      
      cv.min <- cv.return$cv.min
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      lambda <- cv.return$lambda
      basis <- cv.return$basis
    
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
                         restarts=restarts)
      cv.min <- cv.return$cv.min
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
			lambda <- cv.return$lambda
      basis <- cv.return$basis
    
    }
    
  }
  
  est <- crs.default(xz=xz,
                     y=y,
                     degree=degree,
                     segments=segments,
                     include=include,
                     kernel=kernel,
                     lambda=lambda,
                     kernel.type=kernel.type,
                     complexity=complexity,
                     knots=knots,
                     basis=basis,
                     deriv=deriv,
                     data.return=data.return,
                     prune=prune,
                     ...)
  
  est$call <- match.call()
  est$formula <- formula
  est$terms <- mt
  est$xlevels <- .getXlevels(mt, mf)
  est$xz <- xz
  est$y <- y
  est$prune <- prune
  est$cv.min <- cv.min
  est$restarts <- restarts
  
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
    rm(xztmp)

    ## Get evaluation data (newdata) and divide into factors and
    ## numeric.

    Terms <- delete.response(terms(object))
    newdata <- model.frame(Terms,newdata,xlev=object$xlevels)

    ## May 20 - this could be a solution of sorts... the issue is that xz does not have information
    #newdata <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]

    if(!object$kernel) {
      xztmp <- splitFrame(data.frame(newdata))
    } else {
      xztmp <- splitFrame(data.frame(newdata),factor.to.numeric=TRUE)
    }
    xeval <- xztmp$x
    zeval <- xztmp$z
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

      tmp <- predict.factor.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=K,
                                   I=include,
                                   xeval=xeval,
                                   zeval=zeval,
                                   basis=basis,
                                   knots=knots,
                                   prune=prune,
                                   prune.index=prune.index)$fitted.values

      fitted.values <- tmp[,1]
      lwr <- tmp[,2]
      upr <- tmp[,3]
      rm(tmp)

      if(deriv!=0) {

        deriv.mat <- matrix(NA,nrow=NROW(newdata),ncol=NCOL(newdata))
        deriv.mat.lwr <- deriv.mat
        deriv.mat.upr <- deriv.mat
        l <- 1 ## num.z
        m <- 1 ## num.x
        for(i in 1:ncol(newdata)) {
          if(!is.factor(newdata[,i])) {
            tmp <- deriv.factor.spline(x=x,
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
                                       prune.index=prune.index)
            deriv.mat[,i] <- tmp[,1]
            deriv.mat.lwr[,i] <- tmp[,2]
            deriv.mat.upr[,i] <- tmp[,3]
            rm(tmp) 
            m <- m + 1
          } else {
            zevaltmp <- zeval
            zevaltmp[,l] <- factor(rep(levels(newdata[,i])[1],NROW(newdata)),levels=levels(newdata[,i]))
            zpred <- predict.factor.spline(x=x,
                                           y=y,
                                           z=z,
                                           K=K,
                                           I=include,
                                           xeval=xeval,
                                           zeval=zeval,
                                           knots=knots,
                                           basis=basis,
                                           prune=prune,
                                           prune.index=prune.index)$fitted.values

            zpred.base <- predict.factor.spline(x=x,
                                                y=y,
                                                z=z,
                                                K=K,
                                                I=include,
                                                xeval=xeval,
                                                zeval=zevaltmp,
                                                knots=knots,
                                                basis=basis,
                                                prune=prune,
                                                prune.index=prune.index)$fitted.values

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

      kernel.type <- object$kernel.type

      z <- as.matrix(z)
      zeval <- as.matrix(zeval)

      tmp <- predict.kernel.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=K,
                                   lambda=lambda,
                                   kernel.type=kernel.type,
                                   xeval=xeval,
                                   zeval=zeval,
                                   knots=knots,
                                   basis=basis)$fitted.values

      fitted.values <- tmp[,1]
      lwr <- tmp[,2]
      upr <- tmp[,3]
      rm(tmp)

      if(deriv!=0) {

        deriv.mat <- matrix(NA,nrow=NROW(newdata),ncol=NCOL(newdata))
        deriv.mat.lwr <- deriv.mat
        deriv.mat.upr <- deriv.mat
        l <- 1 ## num.z
        m <- 1 ## num.x
        for(i in 1:ncol(newdata)) {
          if(!is.factor(newdata[,i])) {
            tmp <- deriv.kernel.spline(x=x,
                                       y=y,
                                       z=z,
                                       K=K,
                                       lambda=lambda,
                                       kernel.type=kernel.type,
                                       xeval=xeval,
                                       zeval=zeval,
                                       knots=knots,
                                       basis=basis,
                                       deriv.index=m,
                                       deriv=deriv)
            deriv.mat[,i] <- tmp[,1]
            deriv.mat.lwr[,i] <- tmp[,2]
            deriv.mat.upr[,i] <- tmp[,3]
            rm(tmp) 
            m <- m + 1
          } else {
            zevaltmp <- zeval
            zevaltmp[,l] <- as.numeric(factor(rep(levels(newdata[,i])[1],NROW(newdata)),levels=levels(newdata[,i])))
            zpred <- predict.kernel.spline(x=x,
                                           y=y,
                                           z=z,
                                           K=K,
                                           lambda=lambda,
                                           kernel.type=kernel.type,
                                           xeval=xeval,
                                           zeval=zeval,
                                           knots=knots,
                                           basis=basis)$fitted.values

            zpred.base <- predict.kernel.spline(x=x,
                                                y=y,
                                                z=z,
                                                K=K,
                                                lambda=lambda,
                                                kernel.type=kernel.type,
                                                xeval=xeval,
                                                zeval=zevaltmp,
                                                knots=knots,
                                                basis=basis)$fitted.values

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
    cat("\nIndicator Bases/B-spline Bases Regression Spline\n",sep="")
  } else {
    cat("\nKernel Weighting/B-spline Bases Regression Spline\n",sep="")
  }
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
  cat(paste("\nKnot type: ", format(object$knots), sep=""))
  cat(paste("\nModel complexity proxy: ", format(object$complexity), sep=""))
    for(j in 1:object$num.x)
      cat(paste("\nSpline degree/number of segments for ",format(object$xnames[j]),": ",format(object$degree[j]),"/",format(object$segments[j]),sep=""),sep="")
  if(!is.null(object$include)) for(j in 1:length(object$include))
    cat(paste("\nInclusion indicator for ",format(object$znames[j]),": ",format(object$include[j]),sep=""),sep="")
  if(!is.null(object$lambda)) for(j in 1:length(object$lambda))
    cat(paste("\nBandwidth for ",format(object$znames[j]),": ",format(object$lambda[j]),sep=""),sep="")
  cat(paste("\nBasis type: ",format(object$basis),sep=""))
  if(!object$kernel) cat(paste("\nPruning of final model: ",format(ifelse(object$prune,"TRUE","FALSE")),sep=""))
  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  cat(paste("\nRank of model frame: ", format(object$k), sep=""))  
  cat(paste("\nResidual standard error: ", format(sqrt(sum(object$residuals^2)/object$df.residual),digits=4)," on ", format(object$df.residual)," degrees of freedom",sep=""))
  adjusted.r.squared <- 1-(1-object$r.squared)*(length(object$fitted.values)-1)/object$df.residual
  cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4),",   Adjusted R-squared: ",format(adjusted.r.squared,digits=4), sep=""))
  cat(paste("\nCross-validation score: ", format(object$cv.score,digits=8), sep=""))  

  if(sigtest&!object$kernel) {
    cat("\n\nPredictor significance test:\n")
    crs.sigtest(object)
  }

  cat("\n\n")

}

plot.crs <- function(x,
                     mean=FALSE,
                     deriv=FALSE,
                     ci=FALSE,
                     num.eval=100,
                     caption=list("Residuals vs Fitted",
                       "Normal Q-Q Plot",
                       "Scale-Location",
                       "Cook's Distance"),
                     plot.behavior = c("plot","plot-data","data"),
                     ...) {

  plot.behavior <- match.arg(plot.behavior)

  ## We use object below as x is used for data but plot wants
  ## function(x,..)

  object <- x

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  ## Default - basic residual plots

  if(!mean&!deriv) {

    par(mfrow=c(2,2))

    ## Residuals versus fitted

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Working...",console = console)

    plot(fitted(object),
         residuals(object),
         xlab="Fitted Values",
         ylab="Residuals",
         main=caption[[1]],
         ...)

    ## QQ plot

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Working...",console = console)

    std.res <- residuals(object)/sqrt(mean(residuals(object)^2))

    qqnorm(std.res,
           main=caption[[2]])

    qqline(std.res)

    ## Standardized versus fitted

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Working...",console = console)

    plot(fitted(object),
         sqrt(abs(residuals(object,"pearson"))),
         xlab="Fitted Values",
         ylab=as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Standardized Residuals")))),
         main=caption[[3]],
         ...)

    ## Cook's distance \frac{\hat\epsilon_t^2
    ## h_{tt}}{\hat\sigma^2(1-h_{tt})^2k}. Note that this is not
    ## computed for kernel method (what is df.residual when there are
    ## multiple models, etc.?)


    if(!object$kernel) {
      
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush("Working...",console = console)

      sigmasq <- sum(residuals(object)^2)/object$df.residual
      cook <- (residuals(object)^2*object$hatvalues)/(sigmasq*(1-object$hatvalues)^2*object$k)
      plot(cook,
           type = "h",
           main=caption[[4]],
           xlab = "Observation Number",
           ylab = "Cook's Distance",
           ...)
      
    }

    console <- printClear(console)
    console <- printPop(console)
      
  }

  ## Mean

  if(mean) {

    if(!is.null(object$num.z)||(object$num.x>1)) par(mfrow=dim.plot(NCOL(object$xz)))

    mg <- list()

    ## Drawback - data must be first cast outside formula for plot to
    ## work properly (Tristen figured this out so can hunt down issue
    ## later - but this works)

    for(i in 1:NCOL(object$xz)) {

      if(!is.factor(object$xz[,i])) {
        newdata <- matrix(NA,nrow=num.eval,ncol=NCOL(object$xz))
        neval <- num.eval
      } else {
        newdata <- matrix(NA,nrow=length(levels(object$xz[,i])),ncol=NCOL(object$xz))
        neval <- length(levels(object$xz[,i]))
      }

      newdata <- data.frame(newdata)

      if(!is.factor(object$xz[,i])) {
        newdata[,i] <- seq(min(object$xz[,i]),max(object$xz[,i]),length=neval)
      } else {
        newdata[,i] <- factor(levels(object$xz[,i]),levels=levels(object$xz[,i]))
      }

      for(j in (1:NCOL(object$xz))[-i]) {
        if(!is.factor(object$xz[,j])) {
          newdata[,j] <- rep(uocquantile(object$xz[,j],.5),neval)
        } else {
          newdata[,j] <- factor(rep(uocquantile(object$xz[,j],.5),neval),levels=levels(object$xz[,j]))
        }
      }

      newdata <- data.frame(newdata)
      names(newdata) <- names(object$xz)

      if(!ci) {

        ## May 20 - trying to debug plot - issue appears to be that
        ## predict is barfing because formula was used and newdata
        ## does not have similar objects...
        
#    print("Here we are")
#    print(class(newdata[,i]))

        mg[[i]] <- data.frame(newdata[,i],predict(object,newdata=newdata))
        names(mg[[i]]) <- c(names(newdata)[i],"deriv")
        
      } else {
        
        fitted.values <- predict(object,newdata=newdata)
        mg[[i]] <- data.frame(newdata[,i],fitted.values,attr(fitted.values,"lwr"),attr(fitted.values,"upr"))
        names(mg[[i]]) <- c(names(newdata)[i],"mean","lwr","upr")
        
      }
      
      console <- printClear(console)
      console <- printPop(console)

    }
    
    ## Can now add common scale for mean if desired.
    
    if(plot.behavior!="data") {
      
      for(i in 1:NCOL(object$xz)) {
        
        if(!ci) {
          plot(mg[[i]][,1],mg[[i]][,2],
               xlab=names(newdata)[i],
               ylab="Conditional Mean",
               type="l",
               ...)
          
        } else {
          ylim <- c(min(mg[[i]][,-1]),max(mg[[i]][,-1]))
          plot(mg[[i]][,1],mg[[i]][,2],
               xlab=names(newdata)[i],
               ylab="Conditional Mean",
               ylim=ylim,
               type="l",
               ...)
          ## Need to overlay for proper plotting of factor errorbars
          par(new=TRUE)
          plot(mg[[i]][,1],mg[[i]][,3],
               xlab="",
               ylab="",
               ylim=ylim,
               type="l",
               axes=FALSE,
               lty=2,
               col=2,
               ...)
          par(new=TRUE)
          plot(mg[[i]][,1],mg[[i]][,4],
               xlab="",
               ylab="",
               ylim=ylim,
               type="l",
               axes=FALSE,
               lty=2,
               col=2,
               ...)
        }
        
      }

    }
      
    if(plot.behavior!="plot") return(mg)
      
  }
    
  ## deriv

  if(deriv) {

    if(object$deriv > 0) {

      if(!is.null(object$num.z)||(object$num.x>1)) par(mfrow=dim.plot(NCOL(object$xz)))

      rg <- list()

      ## Drawback - data must be first cast outside formula for plot to
      ## work properly (Tristen figured this out so can hunt down issue
      ## later - but this works)

      for(i in 1:NCOL(object$xz)) {

        if(!is.factor(object$xz[,i])) {
          newdata <- matrix(NA,nrow=num.eval,ncol=NCOL(object$xz))
          neval <- num.eval
        } else {
          newdata <- matrix(NA,nrow=length(levels(object$xz[,i])),ncol=NCOL(object$xz))
          neval <- length(levels(object$xz[,i]))
        }

        newdata <- data.frame(newdata)

        if(!is.factor(object$xz[,i])) {
          newdata[,i] <- seq(min(object$xz[,i]),max(object$xz[,i]),length=neval)
        } else {
          newdata[,i] <- factor(levels(object$xz[,i]),levels=levels(object$xz[,i]))
        }

        for(j in (1:NCOL(object$xz))[-i]) {
          if(!is.factor(object$xz[,j])) {
            newdata[,j] <- rep(uocquantile(object$xz[,j],.5),neval)
          } else {
            newdata[,j] <- factor(rep(uocquantile(object$xz[,j],.5),neval),levels=levels(object$xz[,j]))
          }
        }

        newdata <- data.frame(newdata)
        names(newdata) <- names(object$xz)

        if(!ci) {
          
          rg[[i]] <- data.frame(newdata[,i],attr(predict(object,newdata=newdata),"deriv.mat")[,i])
          names(rg[[i]]) <- c(names(newdata)[i],"deriv")
          
        } else {
          
          fitted.values <- predict(object,newdata=newdata)
          rg[[i]] <- data.frame(newdata[,i],
                                attr(predict(object,newdata=newdata),"deriv.mat")[,i],
                                attr(predict(object,newdata=newdata),"deriv.mat.lwr")[,i],
                                attr(predict(object,newdata=newdata),"deriv.mat.upr")[,i])
          names(rg[[i]]) <- c(names(newdata)[i],"deriv","lwr","upr")
          
        }
      
        console <- printClear(console)
        console <- printPop(console)

      }

    } else {

      ## If no deriv given (default=0) issue warning and return

      warning(paste(" derivative plot requested but derivative order is", object$deriv),": specify `deriv=' in crs call",sep="")

    }
    
    ## Can now add common scale for mean if desired.
    
    if(object$deriv > 0) {

      if(plot.behavior!="data") {

        for(i in 1:NCOL(object$xz)) {
          
          if(!ci) {
            plot(rg[[i]][,1],rg[[i]][,2],
                 xlab=names(newdata)[i],
                 ylab=ifelse(!is.factor(newdata[,i]), paste("Order", object$deriv,"Derivative"), "Difference in Levels"),
                 type="l",
                 ...)
            
          } else {
            ylim <- c(min(rg[[i]][,-1]),max(rg[[i]][,-1]))
            plot(rg[[i]][,1],rg[[i]][,2],
                 xlab=names(newdata)[i],
                 ylab=ifelse(!is.factor(newdata[,i]), paste("Order", object$deriv,"Derivative"), "Difference in Levels"),
                 ylim=ylim,
                 type="l",
                 ...)
            ## Need to overlay for proper plotting of factor errorbars
            par(new=TRUE)
            plot(rg[[i]][,1],rg[[i]][,3],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2,
                 ...)
            par(new=TRUE)
            plot(rg[[i]][,1],rg[[i]][,4],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2,
                 ...)
          }
          
        }
        
      }
      
      if(plot.behavior!="plot") return(rg)
      
    }
    
  }

  ## Reset par to 1,1 (can be modified above)
  
  par(mfrow=c(1,1))

}

crs.sigtest <- function(object,...) {

  if(object$kernel) stop(" sigtest is currently available only when kernel=FALSE")

  sg <- list()

  ## Conduct the significance test in order variable by variable

  j.num.x <- 1
  j.num.z <- 1  

  for(i in 1:NCOL(object$xz)) {
    
    if(!is.factor(object$xz[,i])) {
      degree <- object$degree
      degree[j.num.x] <- 0
      model.res <- crs(object$formula,cv="none",degree=degree,include=object$include,basis=object$basis,prune=object$prune,data=eval(object$call$data))
      sg[[i]] <- anova(model.res$model.lm,object$model.lm)
      j.num.x <- j.num.x + 1
    } else {
      include <- object$include
      include[j.num.z] <- 0
      model.res <- crs(object$formula,cv="none",degree=object$degree,include=include,basis=object$basis,prune=object$prune,data=eval(object$call$data))
      sg[[i]] <- anova(model.res$model.lm,object$model.lm)
      j.num.z <- j.num.z + 1      
    }

    cat(paste("Predictor ", format(names(object$xz)[i]), ": Df = ", sg[[i]]$Df[2], ", F = ", format(sg[[i]]$F[2],digits=4), ", Pr(>F) = ", format(sg[[i]][[6]][2],digits=4), "\n", sep=""))
    
  }

}
