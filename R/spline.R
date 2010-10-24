## Generate tensor product spline functions with choice of spline. Can
## be used for training or evaluation data. For continuous datatypes
## uses the B-spline, for factor datatypes uses indicator splines
## obtained from model.matrix(). It also computes derivatives for the
## continuous variables of arbitrary order (issues warning when order
## exceeds degree of spline) with interaction if specified.

prod.spline <- function(x,
                        z=NULL,
                        K=NULL,
                        I=NULL,
                        degree=3,
                        nbreak=2,
                        xeval=NULL,
                        zeval=NULL,
                        basis=c("additive-tensor","additive","tensor","auto"),
                        deriv.index=1,
                        deriv=0) {

  basis <- match.arg(basis)

  if(missing(x) || missing (K)) stop(" must provide x and K")

  ## Care in passing (extra cast) and ensure K is a vector of integers
  ## (K contains the spline degree [integer] for each dimension).

  x <- as.matrix(x)
  K <- round(K) 

  n <- NROW(x)
  num.x <- NCOL(x)
  num.K <- NROW(K)

  if(deriv < 0) stop(" deriv is invalid")
  if(deriv > K[deriv.index]) warning(" deriv order too large, result will be zero")
  if(deriv.index < 1 || deriv.index > num.x) stop(" deriv.index is invalid")

  if(!is.null(z)) {
    z <- data.frame(z)
    num.z <- NCOL(z)
    num.I <- NROW(I)
    if(!is.null(zeval)) {
      zeval <- data.frame(zeval)
    }
  }

  if(is.null(xeval)) {
    xeval <- as.matrix(x)
  } else {
    xeval <- as.matrix(xeval)
    if(NCOL(x)!=NCOL(xeval)) stop(" xeval must be of the same dimension as x")
  }

  if(num.K != num.x) stop(paste(" dimension of x and K incompatible (",num.x,",",num.K,")",sep=""))
  if(!is.null(z) && (num.I != num.z)) stop(paste(" dimension of z and I incompatible (",num.z,",",num.I,")",sep=""))

  if(any(K > 0)||any(I != 0)) {

    tp <- list()

    j <- 1
    for(i in 1:num.x) {
      if(K[i] > 0) {
        if(i==deriv.index) {
          tp[[j]] <- predict(gsl.bs(x[,i,drop=FALSE],degree=K[i],nbreak=nbreak,deriv=deriv,intercept=FALSE),newx=xeval[,i,drop=FALSE])
        } else {
          tp[[j]] <- predict(gsl.bs(x[,i,drop=FALSE],degree=K[i],nbreak=nbreak,intercept=FALSE),newx=xeval[,i,drop=FALSE])
        }
        j <- j+1
      }
    }

    if(!is.null(z)) for(i in 1:num.z) {
      if(I[i] == 1) {
        if(is.null(zeval)) {
          tp[[j]] <- model.matrix(~z[,i])[,-1,drop=FALSE]
        } else {
          tp[[j]] <- model.matrix(~zeval[,i])[,-1,drop=FALSE]
        }
        j <- j+1
      }
    }

    ## When more than one element of K > 0 or I > 0 take all bases
    ## plus tensor product (all interactions), otherwise just the
    ## original bases for the one variable.

    ## [Oct 23 2010, CESG Vancouver] Currently we have additive
    ## splines with (additive) tensor product added via
    ## cbind(P,tensor.prod.model.matrix(tp)) [`functional anova'
    ## setup]. We might simply add a switch at this point to only
    ## return the tensor product as in P<-tensor.prod.model.matrix(tp)
    ## and we are done. Then we have additive, additive with tensor
    ## (`functional anova') and tensor only.

    if(NROW(tp) > 1) {
      ## First create all basis matrices for all continuous predictors
      ## (in essence, additive by default)
      P <- tp[[1]]
      for(i in 2:NROW(tp)) P <- cbind(P,tp[[i]])
      dim.P.no.tensor <- NCOL(P)
      ## Now append tensor to additive if basis==additive-tensor
      if(basis=="additive-tensor") P <- cbind(P,tensor.prod.model.matrix(tp))
      ## Solely tensor if basis==tensor      
      if(basis=="tensor") P <- tensor.prod.model.matrix(tp)
    } else {
      P <- tp[[1]]
      dim.P.no.tensor <- NCOL(P)
    }

  } else {

    ## No relevant continuous or discrete predictors.
    dim.P.no.tensor <- 0
    P <- matrix(rep(1,num.x),num.x,1)

  }

  attr(P,"dim.P.no.tensor") <- dim.P.no.tensor

  return(P)

}

## This function returns the fitted/predicted values for the spline
## regression model with kernel smoothing of the discrete covariates.

predict.kernel.spline <- function(x,
                                  y,
                                  z=NULL,
                                  K,
                                  degree=3,
                                  nbreak=2,
                                  lambda=NULL,
                                  kernel.type=c("nominal","ordinal"),
                                  xeval=NULL,
                                  zeval=NULL,
                                  basis=c("additive-tensor","additive","tensor","auto")){

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")

  basis <- match.arg(basis)
  kernel.type <- match.arg(kernel.type)

  x <- as.matrix(x)

  if(!is.null(z)) z <- as.matrix(z)

  console <- newLineConsole()
  console <- printPush("Working...",console = console)

  if(is.null(z)) {

    ## First no categorical predictor case

    if(any(K > 0)) {

      ## Degree > 0

      P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis)
      model <- lm(y~P)
      if(is.null(xeval)) {
        fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
      } else {
        P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,xeval=xeval,basis=basis)
        fit.spline <- predict(model,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
      }

    } else {

      ## Degree == 0

      model <- lm(y~1)
      if(is.null(xeval)) {
        fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
      } else {
        fit.spline <- predict(model,newdata=data.frame(rep(coef(model),NROW(xeval))),interval="confidence",se.fit=TRUE)
      }
    }

    fit.spline <- cbind(fit.spline[[1]],se=fit.spline[[2]])

    htt <- hatvalues(model)

  } else {

    model <- list()

    ## Categorical predictor case

    n <- NROW(x)

    ## Estimation z information

    z.unique <- uniquecombs(as.matrix(z))
    num.z <- ncol(z.unique)
    ind <-  attr(z.unique,"index")
    ind.vals <-  unique(ind)
    nrow.z.unique <- nrow(z.unique)

    if(any(K > 0)) {

      ## Degree > 0, fitted

      if(is.null(xeval)) {
        fit.spline <- matrix(NA,nrow=n,ncol=4)
        htt <- numeric(length=n)
        for(i in 1:nrow.z.unique) {
          zz <- ind == ind.vals[i]
          L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,kernel.type=kernel.type)
          P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis)
          k <- NCOL(P)
          model.z.unique <- lm(y~P,weights=L)
          model[[i]] <- model.z.unique
          htt[zz] <- hatvalues(model.z.unique)[zz]
          P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,xeval=x[zz,,drop=FALSE],basis=basis)
          tmp <- predict(model.z.unique,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
          fit.spline[zz,] <- cbind(tmp[[1]],se=tmp[[2]])
          rm(tmp)
        }
      } else {

        ## Degree > 0, evaluation

        zeval.unique <- uniquecombs(as.matrix(zeval))
        num.zeval <- ncol(zeval.unique)
        ind.zeval <-  attr(zeval.unique,"index")
        ind.zeval.vals <-  unique(ind.zeval)
        nrow.zeval.unique <- nrow(zeval.unique)

        num.eval <- nrow(zeval)

        fit.spline <- matrix(NA,nrow=num.eval,ncol=4)
        htt <- NULL ## No hatvalues for evaluation
        for(i in 1:nrow.zeval.unique) {
          zz <- ind.zeval == ind.zeval.vals[i]
          L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda,kernel.type=kernel.type)
          P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis)
          k <- NCOL(P)
          model.z.unique <- lm(y~P,weights=L)
          model[[i]] <- model.z.unique
          P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,xeval=xeval[zz,,drop=FALSE],basis=basis)
          tmp <- predict(model.z.unique,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
          fit.spline[zz,] <- cbind(tmp[[1]],se=tmp[[2]])
          rm(tmp)
        }

      }

    } else {

      ## XXX where is predict here?

      ## Degree == 0 (no relevant continuous predictors), train

      if(is.null(xeval)) {

        z.factor <- data.frame(factor(z[,1]))
        k <- ncol(z.factor)
        if(num.z > 1) for(i in 2:num.z) z.factor <- data.frame(z.factor,factor(z[,i]))

        fit.spline <- matrix(NA,nrow=n,ncol=3)
        htt <- numeric(length=n)
        for(i in 1:nrow.z.unique) {
          zz <- ind == ind.vals[i]
          L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,kernel.type=kernel.type)
          model.z.unique <- lm(y~.-1,weights=L,data=data.frame(y,z.factor))
          model[[i]] <- model.z.unique
          htt[zz] <- hatvalues(model.z.unique)[zz]
          fit.spline[zz,] <- fitted(model.z.unique)[zz]
          k <- 0
        }

      } else {

        ## Degree == 0, eval

        z.factor <- data.frame(factor(z[,1]))
        k <- ncol(z.factor)
        if(num.z > 1) for(i in 2:num.z) z.factor <- data.frame(z.factor,factor(z[,i]))

        zeval.unique <- uniquecombs(as.matrix(zeval))
        num.zeval <- ncol(zeval.unique)
        ind.zeval <-  attr(zeval.unique,"index")
        ind.zeval.vals <-  unique(ind.zeval)
        nrow.zeval.unique <- nrow(zeval.unique)

        num.eval <- nrow(zeval)

        fit.spline <- matrix(NA,nrow=num.eval,ncol=3)
        htt <- NULL # no hatvalues for evaluation
        for(i in 1:nrow.zeval.unique) {
          zz <- ind.zeval == ind.zeval.vals[i]
          L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda,kernel.type=kernel.type)
          model.z.unique <- lm(y~.-1,weights=L,data=data.frame(y,z.factor))
          model[[i]] <- model.z.unique
          fit.spline[zz,] <- fitted(model.z.unique)[zz]
          k <- 0
        }

      }

    }

  }

  console <- printClear(console)
  console <- printPop(console)

  return(list(fitted.values=fit.spline,
              df.residual=length(y)-model[[1]]$rank,
              rank=model[[1]]$rank, ## rank same for all models
              model=model,
              hatvalues=htt))

}

## This function returns the gradients of order l and differences in
## levels (order 1 only) for the kernel spline.

deriv.kernel.spline <- function(x,
                                y,
                                z=NULL,
                                K,
                                degree=3,
                                nbreak=2,
                                lambda=NULL,
                                kernel.type=c("nominal","ordinal"),
                                xeval=NULL,
                                zeval=NULL,
                                basis=c("additive-tensor","additive","tensor","auto"),
                                deriv.index=1,
                                deriv=0) {

  if(deriv == 0) stop(" deriv must be greater than zero")

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")

  basis <- match.arg(basis)
  kernel.type <- match.arg(kernel.type)

  x <- as.matrix(x)

  if(!is.null(z)) z <- as.matrix(z)

  if(is.null(z)) {

    ## First no categorical predictor case

    if(K[deriv.index]!=0) {

      P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis)      
      model <- lm(y~P)
      P.deriv <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,xeval=xeval,basis=basis,deriv.index=deriv.index,deriv=deriv)
      dim.P.deriv <- K[deriv.index]
      dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
      dim.P.tensor <- NCOL(P)

      deriv.start <- ifelse(deriv.index!=1,sum(K[1:(deriv.index-1)]),0)+1
      deriv.end <- deriv.start+K[deriv.index]-1

      if(dim.P.tensor > dim.P.deriv+dim.P.no.tensor) {
        deriv.ind.vec <- c(deriv.start:deriv.end, (dim.P.no.tensor+1):dim.P.tensor)
      } else {
        deriv.ind.vec <- 1:dim.P.deriv
      }

      deriv.spline <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]

      vcov.model <- vcov(model)[-1,-1,drop=FALSE]
      se.deriv <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })

    } else {

      deriv.spline <- rep(0,NROW(x))
      se.deriv <- deriv.spline

    }

  } else {

    ## Categorical predictor case

    n <- NROW(x)

    ## Estimation z information

    z.unique <- uniquecombs(as.matrix(z))
    num.z <- ncol(z.unique)
    ind <-  attr(z.unique,"index")
    ind.vals <-  unique(ind)
    nrow.z.unique <- nrow(z.unique)

    if(K[deriv.index]!=0) {

      ## Degree > 0, fitted

      if(is.null(xeval)) {
        deriv.spline <- numeric(length=n)
        se.deriv <- numeric(length=n)        
        for(i in 1:nrow.z.unique) {
          zz <- ind == ind.vals[i]
          L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,kernel.type=kernel.type)
          P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis)          
          k <- NCOL(P)
          model <- lm(y~P,weights=L)
          P.deriv <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,xeval=x[zz,,drop=FALSE],basis=basis,deriv.index=deriv.index,deriv=deriv)
          dim.P.deriv <- K[deriv.index]
          dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
          dim.P.tensor <- NCOL(P)
          deriv.start <- ifelse(deriv.index!=1,sum(K[1:(deriv.index-1)]),0)+1
          deriv.end <- deriv.start+K[deriv.index]-1
          if(dim.P.tensor > dim.P.deriv+dim.P.no.tensor) {
            deriv.ind.vec <- c(deriv.start:deriv.end, (dim.P.no.tensor+1):dim.P.tensor)
          } else {
            deriv.ind.vec <- 1:dim.P.deriv
          }
          deriv.coef <- (coef(model)[-1])[deriv.ind.vec]
          deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%deriv.coef
          vcov.model <- vcov(model)[-1,-1,drop=FALSE]
          se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
        }
      } else {
        ## Evaluation z information

        zeval.unique <- uniquecombs(as.matrix(zeval))
        num.zeval <- ncol(zeval.unique)
        ind.zeval <-  attr(zeval.unique,"index")
        ind.zeval.vals <-  unique(ind.zeval)
        nrow.zeval.unique <- nrow(zeval.unique)

        num.eval <- nrow(zeval)

        deriv.spline <- numeric(length(num.eval))
        se.deriv <- numeric(length=num.eval)        
        for(i in 1:nrow.zeval.unique) {
          zz <- ind.zeval == ind.zeval.vals[i]
          L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda,kernel.type=kernel.type)
          P <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis)
          k <- NCOL(P)
          model <- lm(y~P,weights=L)
          P.deriv <- prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,xeval=xeval[zz,,drop=FALSE],basis=basis,deriv.index=deriv.index,deriv=deriv)
          dim.P.deriv <- K[deriv.index]
          dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
          dim.P.tensor <- NCOL(P)
          deriv.start <- ifelse(deriv.index!=1,sum(K[1:(deriv.index-1)]),0)+1
          deriv.end <- deriv.start+K[deriv.index]-1
          if(dim.P.tensor > dim.P.deriv+dim.P.no.tensor) {
            deriv.ind.vec <- c(deriv.start:deriv.end, (dim.P.no.tensor+1):dim.P.tensor)
          } else {
            deriv.ind.vec <- 1:dim.P.deriv
          }
          deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
          vcov.model <- vcov(model)[-1,-1,drop=FALSE]
          se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
        }
      }

    } else {

      if(is.null(xeval)) {
        deriv.spline <- rep(0,NROW(x))
        se.deriv <- deriv.spline
      } else {
        deriv.spline <- rep(0,NROW(xeval))
        se.deriv <- deriv.spline
      }

    }

  }

  lwr <- deriv.spline - qnorm(0.975)*se.deriv
  upr <- deriv.spline + qnorm(0.975)*se.deriv    

  return(cbind(as.numeric(deriv.spline),lwr, upr))

}

## This function returns the fitted/predicted values using Friedman's
## MARS idea of indicator function bases for categorical variables
## (though Friedman's MARS is much more restrictive than the setup we
## consider here as it uses piecewise linear splines). My additional
## twist is, as for the basis splines, that we allow a variable to not
## enter via a basis of zero length.

predict.factor.spline <- function(x,
                                  y,
                                  z=NULL,
                                  K=NULL,
                                  I=NULL,
                                  degree=3,
                                  nbreak=2,
                                  xeval=NULL,
                                  zeval=NULL,
                                  basis=c("additive-tensor","additive","tensor","auto"),
                                  prune=FALSE,
                                  prune.index=NULL,
                                  trace=0){

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")

  basis <- match.arg(basis)

  ## Cast in case input is not properly cast

  x <- as.matrix(x)
  if(!is.null(xeval)) xeval <- as.matrix(xeval)
  if(!is.null(z)) z <- data.frame(z)
  if(!is.null(zeval)) zeval <- data.frame(zeval)

  console <- newLineConsole()
  console <- printPush("Working...",console = console)

  if(any(K > 0)||any(I>0)) {

    ## Degree > 0

    P <- prod.spline(x=x,z=z,K=K,I=I,degree=degree,nbreak=nbreak,basis=basis)

    if(prune && is.null(prune.index)) {

      ## Pruning via stepwise CV but returning the pruned model only
      ## if the cross-validation score is improved (lower). We create
      ## a data frame so that we can readily determine columns that
      ## have been removed and assign logical values to all columns in
      ## P.
      P.df <- data.frame(P)
      names(P.df) <- paste("P",seq(1,NCOL(P.df)),sep="")
      model <- lm(y~.,data=P.df)
      cv <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      console <- printClear(console)
      console <- printPush("Pruning...",console = console)
      model.pruned <- stepCV(lm(y~.,data=P.df),
                             scope=list(upper=~.,lower=~1),
                             k=log(length(y)),
                             trace=trace)
      cv.pruned <- mean(residuals(model.pruned)^2/(1-hatvalues(model.pruned))^2)
      if(cv.pruned <= cv) {
        IND <- logical()
        for(i in 1:NCOL(P.df)) IND[i] <- any(names(P.df)[i]==names(model.pruned$model[,-1,drop=FALSE]))
        model <- lm(y~P[,IND,drop=FALSE])
      } else {
        warning(" pruned model did not lower cross-validation score, using non-pruned bases")
        IND <- !logical(length=NCOL(P))
        model <- lm(y~P)        
      }
    } else if(prune) {
      ## Pruning, index passed in...
      IND <- prune.index
      model <- lm(y~P[,IND,drop=FALSE])
      cv <- NULL
      cv.pruned <- mean(residuals(model)^2/(1-hatvalues(model))^2)
    } else {
      ## No pruning
      IND <- !logical(length=NCOL(P))
      model <- lm(y~P)
      cv <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      cv.pruned <- NULL
    }      

    if(is.null(xeval)) {
      fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
    } else {
      P <- prod.spline(x=x,z=z,K=K,I=I,degree=degree,nbreak=nbreak,xeval=xeval,zeval=zeval,basis=basis)
      fit.spline <- predict(model,newdata=data.frame(as.matrix(P[,IND,drop=FALSE])),interval="confidence",se.fit=TRUE)
    }

  } else {

    ## Degree == 0, no pruning possible
    IND <- TRUE

    model <- lm(y~1)
    cv <- mean(residuals(model)^2/(1-hatvalues(model))^2) ## Added
    cv.pruned <- NULL
    if(is.null(xeval)) {
      fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
    } else {
      fit.spline <- predict(model,newdata=data.frame(rep(coef(model),NROW(xeval))),interval="confidence",se.fit=TRUE)
    }

  }

  fit.spline <- cbind(fit.spline[[1]],se=fit.spline[[2]])

  console <- printClear(console)
  console <- printPop(console)

  return(list(fitted.values=fit.spline,
              df.residual=model$df.residual,
              rank=model$rank,
              model=model,
              hatvalues=hatvalues(model),
              cv=cv,
              cv.pruned=cv.pruned,
              prune=prune,
              prune.index=IND))

}

## This function returns the fitted/predicted values using Friedman's
## MARS idea of indicator function bases for categorical variables
## (though Friedman's MARS is much more restrictive than the setup we
## consider here as it uses piecewise linear splines). My additional
## twist is, as for the basis splines, that we allow a variable to not
## enter via a basis of zero length.

deriv.factor.spline <- function(x,
                                y,
                                z,
                                K=NULL,
                                I=NULL,
                                degree=3,
                                nbreak=2,
                                xeval=NULL,
                                zeval=NULL,
                                basis=c("additive-tensor","additive","tensor","auto"),
                                deriv.index=1,
                                deriv=0,
                                prune.index=NULL) {

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(deriv == 0) stop(" derivative must be a positive integer")

  basis <- match.arg(basis)

  x <- as.matrix(x)

  if(K[deriv.index]!=0) {

    ## Degree > 0

    ## Estimate model on training data.
    
    P <- prod.spline(x=x,z=z,K=K,I=I,degree=degree,nbreak=nbreak,basis=basis)    
    if(is.null(prune.index)) prune.index <- !logical(NCOL(P))
    model <- lm(y~P[,prune.index,drop=FALSE])

    ## Pad the following for proper handling of pruning

    coef.vec.model <- numeric(length=NCOL(P))
    vcov.mat.model <- matrix(0,nrow=NCOL(P),ncol=NCOL(P))

    coef.vec.model[prune.index] <- coef(model)[-1]
    vcov.mat.model[prune.index,prune.index] <- vcov(model)[-1,-1,drop=FALSE]

    P.deriv <- prod.spline(x=x,z=z,K=K,I=I,degree=degree,nbreak=nbreak,xeval=xeval,zeval=zeval,basis=basis,deriv.index=deriv.index,deriv=deriv)

    dim.P.deriv <- K[deriv.index]
    dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
    dim.P.tensor <- NCOL(P)

    deriv.ind.vec <- logical(length=NCOL(P)) ## All false
    
    deriv.start <- ifelse(deriv.index!=1,sum(K[1:(deriv.index-1)]),0)+1
    deriv.end <- deriv.start+K[deriv.index]-1

    if(dim.P.tensor > dim.P.deriv+dim.P.no.tensor) {
      deriv.ind.vec[c(deriv.start:deriv.end, (dim.P.no.tensor+1):dim.P.tensor)] <- TRUE
    } else {
      deriv.ind.vec[deriv.start:deriv.end] <- TRUE
    }

    deriv.ind.vec <- ifelse(prune.index,deriv.ind.vec,FALSE)

    deriv.spline <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%coef.vec.model[deriv.ind.vec]
    se.deriv <- sapply(1:NROW(P.deriv[,deriv.ind.vec,drop=FALSE]), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.mat.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
    lwr <- deriv.spline - qnorm(0.975)*se.deriv
    upr <- deriv.spline + qnorm(0.975)*se.deriv    

  } else {

    ## Degree == 0

    deriv.spline <- rep(0,NROW(xeval))
    lwr <- deriv.spline
    upr <- deriv.spline

  }

  return(cbind(as.numeric(deriv.spline),lwr, upr))

}

## We use the Sherman-Morrison-Woodbury decomposition to efficiently
## calculate the leave-one-out cross-validation function for
## categorical kernel splines.

cv.kernel.spline <- function(x,
                             y,
                             z=NULL,
                             K,
                             lambda=NULL,
                             z.unique,
                             degree=3,
                             nbreak=2,
                             ind,
                             ind.vals,
                             nrow.z.unique,
                             kernel.type=c("nominal","ordinal"),
                             basis=c("additive-tensor","additive","tensor","auto"),
                             cv.norm=c("L2","L1")) {

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")

  basis <- match.arg(basis)
  kernel.type <- match.arg(kernel.type)
  cv.norm <- match.arg(cv.norm)

  if(is.null(z)) {
    ## No categorical predictors
    if(any(K > 0)) {
      model <- lm(y~prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis))
    } else {
      model <- lm(y~1)
    }
    htt <- hatvalues(model)
    htt <- ifelse(htt == 1, 1-.Machine$double.eps, htt)
    epsilon <- residuals(model)
    cv <- mean(epsilon^2/(1-htt)^2)
  } else {
    ## Categorical predictors
    z <- as.matrix(z)
    num.z <- NCOL(z)
    n <- NROW(y)
    epsilon <- numeric(length=n)
    htt <- numeric(length=n)
    if(any(K > 0)) {
      for(i in 1:nrow.z.unique) {
        zz <- ind == ind.vals[i]
        L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,kernel.type=kernel.type)
        model <- lm(y~prod.spline(x=x,K=K,degree=degree,nbreak=nbreak,basis=basis),weights=L)
        epsilon[zz] <- residuals(model)[zz]
        htt[zz] <- hatvalues(model)[zz]
      }
    } else {
      z.factor <- data.frame(factor(z[,1]))
      if(num.z > 1) for(i in 2:num.z) z.factor <- data.frame(z.factor,factor(z[,i]))
      for(i in 1:nrow.z.unique) {
        zz <- ind == ind.vals[i]
        L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,kernel.type=kernel.type)
        model <- lm(y~.,weights=L,data=data.frame(y,z.factor))
        epsilon[zz] <- residuals(model)[zz]
        htt[zz] <- hatvalues(model)[zz]
      }
    }
    htt <- ifelse(htt == 1, 1-.Machine$double.eps, htt)
    cv <- ifelse(cv.norm=="L2",mean(epsilon^2/(1-htt)^2),mean(abs(epsilon)/abs(1-htt)))
  }

  return(cv)

}

## We use the Sherman-Morrison-Woodbury decomposition to efficiently
## calculate the leave-one-out cross-validation function for factor
## splines.

cv.factor.spline <- function(x,
                             y,
                             z=NULL,
                             K,
                             I=NULL,
                             degree=3,
                             nbreak=2,
                             kernel.type=c("nominal","ordinal"),
                             basis=c("additive-tensor","additive","tensor","auto"),
                             cv.norm=c("L2","L1")) {

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")

  basis <- match.arg(basis)
  kernel.type <- match.arg(kernel.type)
  cv.norm <- match.arg(cv.norm)

  if(!is.null(z)) {
    if(any(K > 0)||any(I > 0)) {
      model <- lm(y~prod.spline(x=x,z=z,K=K,I=I,degree=degree,nbreak=nbreak,basis=basis))
    } else {
      model <- lm(y~1)
    }
  } else {
    if(any(K > 0)) {
      model <- lm(y~prod.spline(x=x,z=z,K=K,I=I,degree=degree,nbreak=nbreak,basis=basis))
    } else {
      model <- lm(y~1)
    }
  }
  htt <- hatvalues(model)
  htt <- ifelse(htt == 1, 1-.Machine$double.eps, htt)
  epsilon <- residuals(model)
  cv <- ifelse(cv.norm=="L2",mean(epsilon^2/(1-htt)^2),mean(abs(epsilon)/abs(1-htt)))

  return(cv)

}
