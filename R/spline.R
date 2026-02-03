## We support multivariate numeric and categorical
## predictors. Multivariate splines are additive and tensor. Can be
## used for training or evaluation data. For continuous datatypes uses
## the B-spline, for factor datatypes uses indicator splines obtained
## from model.matrix(). It also computes derivatives for the
## continuous variables of arbitrary order (issues warning when order
## exceeds degree of spline) with interaction if specified.

## Complexity can be modified via the number of knots (segments) or the
## spline degree (degree)

## Helper to compute hat values from .lm.fit output
hat.from.lm.fit <- function(obj) {
  if(!is.null(obj$qr) && !is.null(obj$qraux) && !is.null(obj$rank)) {
    res <- try(.Call("crs_hat_diag", obj$qr, obj$qraux, as.integer(obj$rank),
                     PACKAGE="crs"), silent=TRUE)
    if(!inherits(res, "try-error")) return(res)
  }
  qr_obj <- list(qr=obj$qr, qraux=obj$qraux, pivot=obj$pivot, tol=obj$tol, rank=obj$rank)
  class(qr_obj) <- "qr"
  hat(qr_obj)
}

prod.spline <- function(x,
                        z=NULL,
                        K=NULL,
                        I=NULL,
                        xeval=NULL,
                        zeval=NULL,
                        knots=c("quantiles","uniform"),
                        basis=c("additive","tensor","glp"),
                        deriv.index=1,
                        deriv=0,
                        ...,
                        display.warnings=TRUE,
                        na.rm) {
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  
  if(missing(x) || missing (K)) stop(" must provide x and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")
  
  ## Additive and glp models have intercept=FALSE in gsl.bs but
  ## intercept=TRUE in lm()
  
  gsl.intercept <- ifelse(basis=="additive" || basis=="glp", FALSE, TRUE)
  
  ## Care in passing (extra cast) and ensure K is a matrix of integers
  ## (K contains the spline degree [integer] for each dimension in
  ## column 1 and segments-1 for each dimension in column 2).
  
  x <- as.matrix(x)
  K <- round(K)
  
  n <- NROW(x)
  num.x <- NCOL(x)
  num.K <- nrow(K)
  
  if(deriv < 0) stop(" deriv is invalid")
  if(deriv > K[deriv.index,1]) if(display.warnings) warning(" deriv order too large, result will be zero")
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
  
  if(any(K[,1] > 0)||any(I != 0)) {
    tp <- list()
    j <- 1
    for(i in 1:num.x) {
      if(K[i,1] > 0) {
        ## nbreak is K[i,2]+1
        if(knots=="uniform") {
          knots.vec <- NULL
        } else {
          ## quantile knots
          knots.vec <- as.numeric(quantile(x[,i,drop=FALSE],probs=seq(0,1,length=(K[i,2]+1))))
          #          if(length(unique(sort(knots.vec))) < length(knots.vec)) {
          ## Correct issue of repeated knots points caused by point
          ## mass data (e.g. knots will be c(0,0,0,1,5), repeated
          ## knots will throw off gsl.bs). This adds a trivial
          ## amount to each knot and is only needed by
          ## gsl.bs(). Otherwise we retain only the unique points
          ## but then the dimension of the spline changes which can
          ## throw off predict etc. Note - there is something odd
          ## about what is produced by quantile as unique does not
          ## work as expected. 1e-20 is too small, 1e-10 works.
          knots.vec <- knots.vec + seq(0,1e-10*(max(x[,i,drop=FALSE])-min(x[,i,drop=FALSE])),length=length(knots.vec))
          #          }
        }
        if((i==deriv.index)&&(deriv!=0)) {
          tp[[j]] <- predict(gsl.bs(x[,i,drop=FALSE],degree=K[i,1],nbreak=(K[i,2]+1),knots=knots.vec,deriv=deriv,intercept=gsl.intercept),newx=xeval[,i,drop=FALSE])
        } else {
          tp[[j]] <- predict(gsl.bs(x[,i,drop=FALSE],degree=K[i,1],nbreak=(K[i,2]+1),knots=knots.vec,intercept=gsl.intercept),newx=xeval[,i,drop=FALSE])
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
    ## When more than one element of K[,1] > 0 or I > 0 take all bases
    ## plus tensor product (all interactions), otherwise just the
    ## original bases for the one variable.
    if(NROW(tp) > 1) {
      ## First create all basis matrices for all continuous predictors
      ## (in essence, additive by default)
      P <- tp[[1]]
      for(i in 2:NROW(tp)) P <- cbind(P,tp[[i]])
      dim.P.no.tensor <- NCOL(P)
      ## Solely tensor if basis==tensor
      if(basis=="tensor") P <- tensor.prod.model.matrix(tp)
      if(basis=="glp") {
        P <- glp.model.matrix(tp)
        if(deriv!=0) {
          P.deriv <- list()
          for(i in 1:length(tp)) P.deriv[[i]] <- matrix(0,1,ncol(tp[[i]]))
          deriv.index <- deriv.index - length(which((K[1:deriv.index,1]==0)))
          while(deriv.index<=0) deriv.index <- deriv.index + 1
          P.deriv[[deriv.index]] <- matrix(NA,1,ncol(tp[[deriv.index]]))
          P[,!is.na(as.numeric(glp.model.matrix(P.deriv)))] <- 0
        }
      }
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

predictKernelSpline <- function(x,
                                y,
                                z=NULL,
                                K,
                                lambda=NULL,
                                is.ordered.z=NULL,
                                xeval=NULL,
                                zeval=NULL,
                                knots=c("quantiles","uniform"),
                                basis=c("additive","tensor","glp"),
                                model.return=FALSE,
                                tau=NULL,
                                weights=NULL,
                                display.warnings=TRUE,
                                display.nomad.progress=TRUE,
                                ...){
  
  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")
  
  if(!is.null(tau)) if(tau <= 0) stop(" tau must be > 0")
  if(!is.null(tau)) if(tau >= 1) stop(" tau must be < 1")
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  if(is.null(is.ordered.z)) stop(" is.ordered.z must be provided")
  
  x <- as.matrix(x)
  
  if(!is.null(z)) z <- as.matrix(z)
  
  console <- newLineConsole()
  if(display.nomad.progress) console <- printPush("Working...",console = console)
  
  model <- NULL ## Returned if model=FALSE and there exist categorical
  ## predictors
  
  if(is.null(z)) {
    
    ## First no categorical predictor case, never reached when called by crs()
    
    if(any(K[,1] > 0)) {
      
      ## Degree > 0
      
      P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
      
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))
          model <- lm(y~P,weights=weights)
        else
          suppressWarnings(model <- rq(y~P,tau=tau,method="fn",weights=weights))
      } else {
        if(is.null(tau))
          model <- lm(y~P-1,weights=weights)
        else
          suppressWarnings(model <- rq(y~P-1,tau=tau,method="fn",weights=weights))
      }
      if(is.null(xeval)) {
        fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
      } else {
        P <- prod.spline(x=x,K=K,xeval=xeval,knots=knots,basis=basis,display.warnings=display.warnings)
        fit.spline <- predict(model,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
      }
      
    } else {
      
      ## Degree == 0
      
      if(is.null(tau))
        model <- lm(y~1,weights=weights)
      else
        suppressWarnings(model <- rq(y~1,tau=tau,method="fn",weights=weights))
      if(is.null(xeval)) {
        fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
      } else {
        fit.spline <- predict(model,newdata=data.frame(rep(coef(model),NROW(xeval))),interval="confidence",se.fit=TRUE)
      }
    }
    
    if(is.null(tau))
      fit.spline <- cbind(fit.spline[[1]],se=fit.spline[[2]])
    else
      fit.spline <- cbind(fit.spline,se=ifelse(NCOL(fit.spline)>1,(fit.spline[,3]-fit.spline[,1])/qnorm(0.975),NA))
    
    if(is.null(tau))
      htt <- hatvalues(model)
    else
      htt <- hat(model$qr)
    
    htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
    
    if(is.null(tau))
      rank <- model$rank
    else
      rank <- NCOL(model$x)
    
  } else {
    
    if(model.return) model <- list()
    
    ## Categorical predictor case
    
    n <- NROW(x)
    
    ## Estimation z information
    
    z.unique <- uniquecombs(as.matrix(z))
    num.z <- ncol(z.unique)
    ind <-  attr(z.unique,"index")
    ind.vals <-  unique(ind)
    nrow.z.unique <- nrow(z.unique)
    
    if(any(K[,1] > 0)) {
      
      ## Degree > 0, fitted
      
      if(is.null(xeval)) {
        fit.spline <- matrix(NA,nrow=n,ncol=4)
        htt <- numeric(length=n)
        P.hat <- numeric(length=n)
        for(i in 1:nrow.z.unique) {
          zz <- ind == ind.vals[i]
          L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
          if(!is.null(weights)) L <- weights*L
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
          k <- NCOL(P)
          if(basis=="additive" || basis=="glp") {
            if(is.null(tau))
              model.z.unique <- lm(y~P,weights=L)
            else
              suppressWarnings(model.z.unique <- rq(y~P,weights=L,tau=tau,method="fn"))
            model.z.unique.hat <- lm(y~P,weights=L)
          } else {
            if(is.null(tau))
              model.z.unique <- lm(y~P-1,weights=L)
            else
              suppressWarnings(model.z.unique <- rq(y~P-1,weights=L,tau=tau,method="fn"))
            model.z.unique.hat <- lm(y~P-1,weights=L)
          }
          if(model.return) model[[i]] <- model.z.unique
          if(is.null(tau))
            htt[zz] <- hatvalues(model.z.unique)[zz]
          else
            htt[zz] <- hatvalues(model.z.unique.hat)[zz]
          
          P.hat[zz] <- sum(L)
          P <- prod.spline(x=x,K=K,xeval=x[zz,,drop=FALSE],knots=knots,basis=basis,display.warnings=display.warnings)
          tmp <- predict(model.z.unique,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
          
          if(is.null(tau))
            fit.spline[zz,] <- cbind(tmp[[1]],tmp[[2]])
          else
            fit.spline[zz,] <- cbind(tmp,(tmp[,3]-tmp[,1])/qnorm(0.975))
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
        P.hat <- NULL
        for(i in 1:nrow.zeval.unique) {
          zz <- ind.zeval == ind.zeval.vals[i]
          L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
          if(!is.null(weights)) L <- weights*L
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
          k <- NCOL(P)
          if(basis=="additive" || basis=="glp") {
            if(is.null(tau))
              model.z.unique <- lm(y~P,weights=L)
            else
              suppressWarnings(model.z.unique <- rq(y~P,weights=L,tau=tau,method="fn"))
          } else {
            if(is.null(tau))
              model.z.unique <- lm(y~P-1,weights=L)
            else
              suppressWarnings(model.z.unique <- rq(y~P-1,weights=L,tau=tau,method="fn"))
          }
          if(model.return) model[[i]] <- model.z.unique
          P <- prod.spline(x=x,K=K,xeval=xeval[zz,,drop=FALSE],knots=knots,basis=basis,display.warnings=display.warnings)
          tmp <- predict(model.z.unique,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
          
          if(is.null(tau))
            fit.spline[zz,] <- cbind(tmp[[1]],tmp[[2]])
          else
            fit.spline[zz,] <- cbind(tmp,(tmp[,3]-tmp[,1])/qnorm(0.975))
          
          rm(tmp)
        }
        
      }
    } else {
      
      ## Degree == 0 (no relevant continuous predictors), train
      
      if(is.null(xeval)) {
        fit.spline <- matrix(NA,nrow=n,ncol=4)
        htt <- numeric(length=n)
        P.hat <- numeric(length=n)
        x.intercept <- rep(1,n)
        for(i in 1:nrow.z.unique) {
          zz <- ind == ind.vals[i]
          L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
          if(!is.null(weights)) L <- weights*L
          k <- 0
          ## Whether we use additive, glp, or tensor products, this
          ## model has no continuous predictors hence the intercept is
          ## the parameter that may shift with the categorical
          ## predictors
          if(is.null(tau))
            model.z.unique <- lm(y~x.intercept-1,weights=L)
          else
            suppressWarnings(model.z.unique <- rq(y~x.intercept-1,weights=L,tau=tau,method="fn"))
          model.z.unique.hat <- lm(y~x.intercept-1,weights=L)
          if(model.return) model[[i]] <- model.z.unique
          if(is.null(tau))
            htt[zz] <- hatvalues(model.z.unique)[zz]
          else
            htt[zz] <- hatvalues(model.z.unique.hat)[zz]
          P.hat[zz] <- sum(L)
          tmp <- predict(model.z.unique,newdata=data.frame(x.intercept=x.intercept[zz]),interval="confidence",se.fit=TRUE)
          
          if(is.null(tau))
            fit.spline[zz,] <- cbind(tmp[[1]],tmp[[2]])
          else
            fit.spline[zz,] <- cbind(tmp,(tmp[,3]-tmp[,1])/qnorm(0.975))
          
          rm(tmp)
        }
      } else {
        
        ## Degree == 0 (no relevant continuous predictors), evaluation
        
        zeval.unique <- uniquecombs(as.matrix(zeval))
        num.zeval <- ncol(zeval.unique)
        ind.zeval <-  attr(zeval.unique,"index")
        ind.zeval.vals <-  unique(ind.zeval)
        nrow.zeval.unique <- nrow(zeval.unique)
        
        num.eval <- nrow(zeval)
        
        fit.spline <- matrix(NA,nrow=num.eval,ncol=4)
        htt <- NULL ## No hatvalues for evaluation
        P.hat <- NULL
        x.intercept <- rep(1,n)
        for(i in 1:nrow.zeval.unique) {
          zz <- ind.zeval == ind.zeval.vals[i]
          L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
          if(!is.null(weights)) L <- weights*L
          k <- 0
          if(is.null(tau))
            model.z.unique <- lm(y~x.intercept-1,weights=L)
          else
            suppressWarnings(model.z.unique <- rq(y~x.intercept-1,weights=L,tau=tau,method="fn"))
          if(model.return) model[[i]] <- model.z.unique
          tmp <- predict(model.z.unique,newdata=data.frame(x.intercept=rep(1,num.eval)[zz]),interval="confidence",se.fit=TRUE)
          
          if(is.null(tau))
            fit.spline[zz,] <- cbind(tmp[[1]],tmp[[2]])
          else
            fit.spline[zz,] <- cbind(tmp,(tmp[,3]-tmp[,1])/qnorm(0.975))
          
          rm(tmp)
        }
      }
      
    }
    
    if(is.null(tau))
      rank <- model.z.unique$rank ## same for all models
    else
      rank <- NCOL(model.z.unique$x) ## same for all models
  }
  
  console <- printClear(console)
  console <- printPop(console)
  
  ## Need to return kernel probability estimates. The kernel function
  ## we use does not sum to one so the probability estimates will not
  ## be proper (will not sum to one), so we simply renormalize by the
  ## sum of the unique probabilities. However, when lambda=1 for all
  ## categorical predictors there is only one unique probability value
  ## and the non-proper probability estimates will all equal one so we
  ## trap this case.
  
  P.hat <- P.hat/(sum(unique(P.hat/n))*n)
  P.hat <- ifelse(P.hat==1,1/nrow.z.unique,P.hat)
  
  return(list(fitted.values=fit.spline,
              df.residual=length(y)-rank,
              rank=rank,
              model=model,
              hatvalues=htt,
              P.hat=P.hat,
              tau=tau))
  
}

## This function returns the gradients of order l and differences in
## levels (order 1 only) for the kernel spline.

derivKernelSpline <- function(x,
                              y,
                              z=NULL,
                              K,
                              lambda=NULL,
                              is.ordered.z=NULL,
                              xeval=NULL,
                              zeval=NULL,
                              knots=c("quantiles","uniform"),
                              basis=c("additive","tensor","glp"),
                              deriv.index=1,
                              deriv=0,
                              tau=NULL,
                              weights=NULL,
                              display.warnings=TRUE,
                              ...) {
  
  if(deriv == 0) stop(" deriv must be greater than zero")
  
  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  if(is.null(is.ordered.z)) stop(" is.ordered.z must be provided")
  
  x <- as.matrix(x)
  
  ## Univariate additive spline bases have one less column than
  ## univariate tensor spline bases. This is used only for setting
  ## appropriate columns for derivative computation. We also need to
  ## set the segments to 0 when the degree is zero, again only for
  ## derivative computation when using an additive basis.
  
  if(basis=="additive" || basis=="glp") {
    K.additive <- K
    K.additive[,2] <- ifelse(K[,1]==0,0,K[,2])
    K.additive[,1] <- ifelse(K[,1]>0,K[,1]-1,K[,1])
  }
  
  if(!is.null(z)) z <- as.matrix(z)
  
  if(is.null(z)) {
    
    ## First no categorical predictor case (never reached by crs)
    
    if(K[deriv.index,1]!=0) {
      
      P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
      P.deriv <- prod.spline(x=x,K=K,xeval=xeval,knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv,display.warnings=display.warnings)
      dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
      dim.P.tensor <- NCOL(P)
      
      if(basis=="additive") {
        if(is.null(tau))
          model <- lm(y~P,weights=weights)
        else
          suppressWarnings(model <- rq(y~P,tau=tau,method="fn",weights=weights))
        
        dim.P.deriv <- sum(K.additive[deriv.index,])
        deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)
        deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
        deriv.ind.vec <- max(1,deriv.start:deriv.end - length(which(K[,1]==0)))
        deriv.spline <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
        
        if(is.null(tau))
          vcov.model <- vcov(model)[-1,-1,drop=FALSE]
        else
          suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
        
        se.deriv <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
      } else if(basis=="tensor") {
        if(is.null(tau))
          model <- lm(y~P-1,weights=weights)
        else
          suppressWarnings(model <- rq(y~P-1,tau=tau,method="fn",weights=weights))
        
        deriv.spline <- P.deriv%*%coef(model)
        
        if(is.null(tau))
          vcov.model <- vcov(model)
        else
          suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov)
        
        se.deriv <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
      } else if(basis=="glp") {
        if(is.null(tau))
          model <- lm(y~P,weights=weights)
        else
          suppressWarnings(model <- rq(y~P,tau=tau,method="fn",weights=weights))
        deriv.spline <- P.deriv%*%coef(model)[-1]
        
        if(is.null(tau))
          vcov.model <- vcov(model)[-1,-1,drop=FALSE]
        else
          suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
        
        se.deriv <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
      }
      
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
    
    if(K[deriv.index,1]!=0) {
      
      ## Degree > 0, fitted
      
      if(is.null(xeval)) {
        deriv.spline <- numeric(length=n)
        se.deriv <- numeric(length=n)
        for(i in 1:nrow.z.unique) {
          zz <- ind == ind.vals[i]
          L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
          if(!is.null(weights)) L <- weights*L
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
          P.deriv <- prod.spline(x=x,K=K,xeval=x[zz,,drop=FALSE],knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv,display.warnings=display.warnings)
          k <- NCOL(P)
          dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
          dim.P.tensor <- NCOL(P)
          
          if(basis=="additive") {
            if(is.null(tau))
              model <- lm(y~P,weights=L)
            else
              suppressWarnings(model <- rq(y~P,tau=tau,method="fn",weights=L))
            dim.P.deriv <- sum(K.additive[deriv.index,])
            deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)
            deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
            deriv.ind.vec <- deriv.start:deriv.end
            deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
            if(is.null(tau))
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
            
            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
          } else if(basis=="tensor") {
            if(is.null(tau))
              model <- lm(y~P-1,weights=L)
            else
              suppressWarnings(model <- rq(y~P-1,tau=tau,method="fn",weights=L))
            
            deriv.spline[zz] <- P.deriv%*%coef(model)
            
            if(is.null(tau))
              vcov.model <- vcov(model)
            else
              suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov)
            
            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
          } else if(basis=="glp") {
            if(is.null(tau))
              model <- lm(y~P,weights=L)
            else
              suppressWarnings(model <- rq(y~P,tau=tau,method="fn",weights=L))
            
            deriv.spline[zz] <- P.deriv%*%coef(model)[-1]
            
            if(is.null(tau))
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
            
            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
          }
          
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
          L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
          if(!is.null(weights)) L <- weights*L
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
          P.deriv <- prod.spline(x=x,K=K,xeval=xeval[zz,,drop=FALSE],knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv,display.warnings=display.warnings)
          k <- NCOL(P)
          dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
          dim.P.tensor <- NCOL(P)
          
          if(basis=="additive") {
            if(is.null(tau))
              model <- lm(y~P,weights=L)
            else
              suppressWarnings(model <- rq(y~P,weights=L,tau=tau,method="fn"))
            
            dim.P.deriv <- sum(K.additive[deriv.index,])
            deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)
            deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
            deriv.ind.vec <- deriv.start:deriv.end
            deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
            if(is.null(tau))
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
            
            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
          } else if(basis=="tensor") {
            if(is.null(tau))
              model <- lm(y~P-1,weights=L)
            else
              suppressWarnings(model <- rq(y~P-1,weights=L,tau=tau,method="fn"))
            deriv.spline[zz] <- P.deriv%*%coef(model)
            
            if(is.null(tau))
              vcov.model <- vcov(model)
            else
              suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov)
            
            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
          } else if(basis=="glp") {
            if(is.null(tau))
              model <- lm(y~P,weights=L)
            else
              suppressWarnings(model <- rq(y~P,weights=L,tau=tau,method="fn"))
            deriv.spline[zz] <- P.deriv%*%coef(model)[-1]
            
            if(is.null(tau))
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              suppressWarnings(vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
            
            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
          }
          
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
## consider here as it uses piece-wise linear splines). My additional
## twist is, as for the basis splines, that we allow a variable to not
## enter via a basis of zero length.

preditFactorSpline <- function(x,
                               y,
                               z=NULL,
                               K=NULL,
                               I=NULL,
                               xeval=NULL,
                               zeval=NULL,
                               knots=c("quantiles","uniform"),
                               basis=c("additive","tensor","glp"),
                               prune=FALSE,
                               prune.index=NULL,
                               trace=0,
                               tau=NULL,
                               weights=NULL,
                               display.warnings=TRUE,
                               display.nomad.progress=TRUE,
                               ...){
  
  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")
  
  if(!is.null(tau)) if(tau <= 0) stop(" tau must be > 0")
  if(!is.null(tau)) if(tau >= 1) stop(" tau must be < 1")
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  
  ## Cast in case input is not properly cast
  
  x <- as.matrix(x)
  if(!is.null(xeval)) xeval <- as.matrix(xeval)
  if(!is.null(z)) z <- data.frame(z)
  if(!is.null(zeval)) zeval <- data.frame(zeval)
  
  console <- newLineConsole()
  if(display.nomad.progress) console <- printPush("Working...",console = console)
  
  if(any(K[,1] > 0)||any(I>0)) {
    
    ## Degree > 0
    
    P <- prod.spline(x=x,z=z,K=K,I=I,knots=knots,basis=basis,display.warnings=display.warnings)
    
    if(prune && is.null(prune.index)) {
      
      ## Pruning via step-wise CV but returning the pruned model only
      ## if the cross-validation score is improved (lower). We create
      ## a data frame so that we can readily determine columns that
      ## have been removed and assign logical values to all columns in
      ## P.
      
      ## Note - this code is not reachable by crs() since pruning and
      ## regression quantiles is not supported by stepCV (currently we
      ## test and stop())
      
      P.df <- data.frame(P)
      names(P.df) <- paste("P",seq(1,NCOL(P.df)),sep="")
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))
          model <- lm(y~.,data=P.df,weights=weights)
        else
          suppressWarnings(model <- rq(y~.,data=P.df,tau=tau,method="fn",weights=weights))
      } else {
        if(is.null(tau))
          model <- lm(y~.-1,data=P.df,weights=weights)
        else
          suppressWarnings(model <- rq(y~.-1,data=P.df,tau=tau,method="fn",weights=weights))
      }
      if(is.null(tau))
        cv <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      else
        suppressWarnings(cv <- cv.rq(model,tau=tau,weights=weights))
      console <- printClear(console)
      if(display.nomad.progress) console <- printPush("Pruning...",console = console)
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))
          model.pruned <- stepCV(lm(y~.,data=P.df,weights=weights),
                                 scope=list(upper=~.,lower=~1),
                                 k=log(length(y)),
                                 trace=trace,
                                 display.warnings=display.warnings)
        else
          suppressWarnings(model.pruned <- stepCV(rq(y~.,data=P.df,tau=tau,method="fn",weights=weights),
                                                  scope=list(upper=~.,lower=~1),
                                                  k=log(length(y)),
                                                  trace=trace,
                                                  display.warnings=display.warnings))
        
      } else {
        if(is.null(tau))
          model.pruned <- stepCV(lm(y~.-1,data=P.df,weights=weights),
                                 scope=list(upper=~.,lower=~1),
                                 k=log(length(y)),
                                 trace=trace,
                                 display.warnings=display.warnings)
        else
          suppressWarnings(model.pruned <- stepCV(rq(y~.-1,data=P.df,tau=tau,method="fn",weights=weights),
                                                  scope=list(upper=~.,lower=~1),
                                                  k=log(length(y)),
                                                  trace=trace,
                                                  display.warnings=display.warnings))
        
      }
      if(is.null(tau))
        cv.pruned <- mean(residuals(model.pruned)^2/(1-hatvalues(model.pruned))^2)
      else
        suppressWarnings(cv.pruned <- cv.rq(model.pruned,tau=tau,weights=weights))
      
      if(cv.pruned <= cv) {
        IND <- logical()
        for(i in 1:NCOL(P.df)) IND[i] <- any(names(P.df)[i]==names(model.pruned$model[,-1,drop=FALSE]))
        if(basis=="additive" || basis=="glp") {
          if(is.null(tau))
            model <- lm(y~P[,IND,drop=FALSE],weights=weights)
          else
            suppressWarnings(model <- rq(y~P[,IND,drop=FALSE],tau=tau,method="fn",weights=weights))
        } else {
          if(is.null(tau))
            model <- lm(y~P[,IND,drop=FALSE]-1,weights=weights)
          else
            suppressWarnings(model <- rq(y~P[,IND,drop=FALSE]-1,tau=tau,method="fn",weights=weights))
        }
      } else {
        if(display.warnings) warning(" pruned model did not lower cross-validation score, using non-pruned bases")
        IND <- !logical(length=NCOL(P))
        if(basis=="additive" || basis=="glp") {
          if(is.null(tau))
            model <- lm(y~P,weights=weights)
          else
            suppressWarnings(model <- rq(y~P,tau=tau,method="fn",weights=weights))
        } else {
          if(is.null(tau))
            model <- lm(y~P-1,weights=weights)
          else
            suppressWarnings(model <- rq(y~P-1,tau=tau,method="fn",weights=weights))
        }
      }
    } else if(prune) {
      ## Pruning, index passed in...
      IND <- prune.index
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))
          model <- lm(y~P[,IND,drop=FALSE],weights=weights)
        else
          suppressWarnings(model <- rq(y~P[,IND,drop=FALSE],tau=tau,method="fn",weights=weights))
      } else {
        if(is.null(tau))
          model <- lm(y~P[,IND,drop=FALSE]-1,weights=weights)
        else
          suppressWarnings(model <- rq(y~P[,IND,drop=FALSE]-1,tau=tau,method="fn",weights=weights))
      }
      cv <- NULL
      if(is.null(tau))
        cv.pruned <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      else
        suppressWarnings(cv.pruned <- cv.rq(model,tau=tau,weights=weights))
    } else {
      ## No pruning, default case
      IND <- !logical(length=NCOL(P))
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))
          model <- lm(y~P,weights=weights)
        else
          suppressWarnings(model <- rq(y~P,tau=tau,method="fn",weights=weights))
      } else {
        if(is.null(tau))
          model <- lm(y~P-1,weights=weights)
        else
          suppressWarnings(model <- rq(y~P-1,tau=tau,method="fn",weights=weights,x=TRUE))
      }
      cv.pruned <- NULL
      if(is.null(tau)) {
        cv <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      } else {
        if(basis=="additive" || basis=="glp")
          model.hat <- lm(y~P,weights=weights)
        else
          model.hat <- lm(y~P-1,weights=weights)
        htt <- hat(model.hat$qr)
        ## Note - this is defined in util.R so if you modify there you must modify here also
        if(is.null(weights))
          cv <- mean(check.function(residuals(model),tau)/(1-htt)^(1/sqrt(tau*(1-tau))))
        else
          cv <- mean(check.function(residuals(model)*sqrt(weights),tau)/(1-htt)^(1/sqrt(tau*(1-tau))))
      }
    }
    
    if(is.null(xeval)) {
      fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
    } else {
      P <- prod.spline(x=x,z=z,K=K,I=I,xeval=xeval,zeval=zeval,knots=knots,basis=basis,display.warnings=display.warnings)
      fit.spline <- predict(model,newdata=data.frame(as.matrix(P[,IND,drop=FALSE])),interval="confidence",se.fit=TRUE)
    }
    
  } else {
    
    ## Degree == 0, no pruning possible
    IND <- TRUE
    
    if(is.null(tau)) {
      model <- lm(y~1,weights=weights)
      cv <- mean(residuals(model)^2/(1-hatvalues(model))^2) ## Added
    } else {
      suppressWarnings(model <- rq(y~1,tau=tau,method="fn",weights=weights))
      model.hat <- lm(y~1,weights=weights)
      htt <- hat(model.hat$qr)
      ## Note - this is defined in util.R so if you modify there you must modify here also
      if(is.null(weights))
        cv <- mean(check.function(residuals(model),tau)/(1-htt)^(1/sqrt(tau*(1-tau))))
      else
        cv <- mean(check.function(residuals(model)*sqrt(weights),tau)/(1-htt)^(1/sqrt(tau*(1-tau))))
    }
    
    cv.pruned <- NULL
    if(is.null(xeval)) {
      fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
    } else {
      fit.spline <- predict(model,newdata=data.frame(rep(coef(model),NROW(xeval))),interval="confidence",se.fit=TRUE)
    }
    
  }
  
  if(is.null(tau))
    fit.spline <- cbind(fit.spline[[1]],se=fit.spline[[2]])
  else
    fit.spline <- cbind(fit.spline,se=ifelse(NCOL(fit.spline)>1,(fit.spline[,3]-fit.spline[,1])/qnorm(0.975),NA))
  
  console <- printClear(console)
  console <- printPop(console)
  
  if(is.null(tau))
    htt <- hatvalues(model)
  else
    htt <- hatvalues(model.hat)
  
  return(list(fitted.values=fit.spline,
              df.residual=model$df.residual,
              rank=model$rank,
              model=model,
              hatvalues=htt,
              cv=cv,
              cv.pruned=cv.pruned,
              prune=prune,
              prune.index=IND,
              tau=tau))
  
}

## This function returns the fitted/predicted values using Friedman's
## MARS idea of indicator function bases for categorical variables
## (though Friedman's MARS is much more restrictive than the setup we
## consider here as it uses piece-wise linear splines). My additional
## twist is, as for the basis splines, that we allow a variable to not
## enter via a basis of zero length.

derivFactorSpline <- function(x,
                              y,
                              z,
                              K=NULL,
                              I=NULL,
                              xeval=NULL,
                              zeval=NULL,
                              knots=c("quantiles","uniform"),
                              basis=c("additive","tensor","glp"),
                              deriv.index=1,
                              deriv=0,
                              prune.index=NULL,
                              tau=NULL,
                              weights=NULL,
                              display.warnings=TRUE,
                              ...) {
  
  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(deriv == 0) stop(" derivative must be a positive integer")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  
  x <- as.matrix(x)
  
  ## Univariate additive spline bases have one less column than
  ## univariate tensor spline bases. This is used only for setting
  ## appropriate columns for derivative computation. We also need to
  ## set the segments to 0 when the degree is zero, again only for
  ## derivative computation when using an additive basis.
  
  if(basis=="additive" || basis=="glp") {
    K.additive <- K
    K.additive[,2] <- ifelse(K[,1]==0,0,K[,2])
    K.additive[,1] <- ifelse(K[,1]>0,K[,1]-1,K[,1])
  }
  if(K[deriv.index,1]!=0) {
    
    ## Degree > 0
    
    ## Estimate model on training data.
    
    P <- prod.spline(x=x,z=z,K=K,I=I,knots=knots,basis=basis,display.warnings=display.warnings)
    P.deriv <- prod.spline(x=x,z=z,K=K,I=I,xeval=xeval,zeval=zeval,knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv,display.warnings=display.warnings)
    
    dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
    dim.P.tensor <- NCOL(P)
    deriv.ind.vec <- logical(length=NCOL(P)) ## All false
    
    if(is.null(prune.index)) prune.index <- !logical(NCOL(P))
    
    ## Pad the following for proper handling of pruning
    
    coef.vec.model <- numeric(length=NCOL(P))
    vcov.mat.model <- matrix(0,nrow=NCOL(P),ncol=NCOL(P))
    if(basis=="additive") {
      if(is.null(tau))
        model <- lm(y~P[,prune.index,drop=FALSE],weights=weights)
      else
        suppressWarnings(model <- rq(y~P[,prune.index,drop=FALSE],tau=tau,weights=weights,method="fn"))
      
      coef.vec.model[prune.index] <- coef(model)[-1]
      
      if(is.null(tau))
        vcov.mat.model[prune.index,prune.index] <- vcov(model)[-1,-1,drop=FALSE]
      else
        suppressWarnings(vcov.mat.model[prune.index,prune.index] <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
      
      dim.P.deriv <- sum(K.additive[deriv.index,])
      deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)
      deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
      deriv.ind.vec[deriv.start:deriv.end] <- TRUE
      deriv.ind.vec <- ifelse(prune.index,deriv.ind.vec,FALSE)
    } else if(basis=="tensor") {
      if(is.null(tau))
        model <- lm(y~P[,prune.index,drop=FALSE]-1,weights=weights)
      else
        suppressWarnings(model <- rq(y~P[,prune.index,drop=FALSE]-1,tau=tau,weights=weights,method="fn"))
      coef.vec.model[prune.index] <- coef(model)
      
      if(is.null(tau))
        vcov.mat.model[prune.index,prune.index] <- vcov(model)
      else
        suppressWarnings(vcov.mat.model[prune.index,prune.index] <- summary(model,covariance=TRUE)$cov)
      
      deriv.ind.vec[1:dim.P.tensor] <- TRUE
      deriv.ind.vec <- ifelse(prune.index,deriv.ind.vec,FALSE)
    } else if(basis=="glp") {
      if(is.null(tau))
        model <- lm(y~P[,prune.index,drop=FALSE],weights=weights)
      else
        suppressWarnings(model <- rq(y~P[,prune.index,drop=FALSE],tau=tau,weights=weights,method="fn"))
      coef.vec.model[prune.index] <- coef(model)[-1]
      
      if(is.null(tau))
        vcov.mat.model[prune.index,prune.index] <- vcov(model)[-1,-1,drop=FALSE]
      else
        suppressWarnings(vcov.mat.model[prune.index,prune.index] <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE])
      
      deriv.ind.vec[1:dim.P.tensor] <- TRUE
      deriv.ind.vec <- ifelse(prune.index,deriv.ind.vec,FALSE)
    }
    
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

## The following function is a wrapper of cv.kernel.spline to
## handle the situation when knots="auto". It will call
## cv.kernel.spline two times and return the minimum with
## the optimal knots.
cv.kernel.spline.wrapper <- function(x,
                                     y,
                                     z=NULL,
                                     K,
                                     lambda=NULL,
                                     z.unique,
                                     ind,
                                     ind.vals,
                                     ind.list=NULL,
                                     nrow.z.unique,
                                     is.ordered.z=NULL,
                                     knots=c("quantiles","uniform","auto"),
                                     basis=c("additive","tensor","glp"),
                                     cv.func=c("cv.ls","cv.gcv","cv.aic"),
                                     cv.df.min=1,
                                     tau=NULL,
                                     weights=NULL,
                                     singular.ok=FALSE,
                                     display.warnings=TRUE) {
  
  knots.opt <- knots;
  
  if(knots == "auto") {
    
    knots.opt <- "quantiles"
    
    cv <- cv.kernel.spline(x=x,
                           y=y,
                           z=z,
                           K=K,
                           lambda=lambda,
                           z.unique=z.unique,
                           ind=ind,
                           ind.vals=ind.vals,
                           ind.list=ind.list,
                           nrow.z.unique=nrow.z.unique,
                           is.ordered.z=is.ordered.z,
                           knots="quantiles",
                           basis=basis,
                           cv.func=cv.func,
                           cv.df.min=cv.df.min,
                           tau=tau,
                           weights=weights,
                           singular.ok=singular.ok,
                           display.warnings=display.warnings)
    
    cv.uniform <- cv.kernel.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=K,
                                   lambda=lambda,
                                   z.unique=z.unique,
                                   ind=ind,
                                   ind.vals=ind.vals,
                                   ind.list=ind.list,
                                   nrow.z.unique=nrow.z.unique,
                                   is.ordered.z=is.ordered.z,
                                   knots="uniform",
                                   basis=basis,
                                   cv.func=cv.func,
                                   cv.df.min=cv.df.min,
                                   tau=tau,
                                   weights=weights,
                                   singular.ok=singular.ok,
                                   display.warnings=display.warnings)
    if(cv > cv.uniform) {
      cv <- cv.uniform
      knots.opt <- "uniform"
    }
    
  } else {
    
    cv <- cv.kernel.spline(x=x,
                           y=y,
                           z=z,
                           K=K,
                           lambda=lambda,
                           z.unique=z.unique,
                           ind=ind,
                           ind.vals=ind.vals,
                           ind.list=ind.list,
                           nrow.z.unique=nrow.z.unique,
                           is.ordered.z=is.ordered.z,
                           knots=knots,
                           basis=basis,
                           cv.func=cv.func,
                           cv.df.min=cv.df.min,
                           tau=tau,
                           weights=weights,
                           singular.ok=singular.ok,
                           display.warnings=display.warnings)
    
  }
  
  attr(cv, "knots.opt") <- knots.opt
  
  return(cv)
  
}
## We use the Sherman-Morrison-Woodbury decomposition to efficiently
## calculate the leave-one-out cross-validation function for
## categorical kernel splines.

## June 24 2011 - replaced lm() with model <- lm.fit/lm.wfit and
## hat(model$qr) both here and in cv.factor.spline. Here it reduces
## runtime by 20-30%. But more importantly lm.fit is more `robust'
## than lsfit (lm.fit is the `workhorse' of lm, lsfit calls LAPACK
## code). Note that it is noticeable as it returns a larger cv value
## for more complicated problems which is naturally desirable.

# cv.kernel.spline <- function(x,
#                              y,
#                              z=NULL,
#                              K,
#                              lambda=NULL,
#                              z.unique,
#                              ind,
#                              ind.vals,
#                              nrow.z.unique,
#                              is.ordered.z=NULL,
#                              knots=c("quantiles","uniform"),
#                              basis=c("additive","tensor","glp"),
#                              cv.func=c("cv.ls","cv.gcv","cv.aic"),
#                              cv.df.min=1,
#                              tau=NULL,
#                              weights=NULL,
#                              singular.ok=FALSE,
#                              display.warnings=TRUE) {
#   
#   if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
#   
#   if(!is.matrix(K)) stop(" K must be a two-column matrix")
#   
#   basis <- match.arg(basis)
#   if(is.null(is.ordered.z)) stop(" is.ordered.z must be provided")
#   knots <- match.arg(knots)
#   cv.func <- match.arg(cv.func)
#   
#   ## Without computing P, compute the number of columns that P would
#   ## be and if degrees of freedom is 1 or less, return a large penalty.
#   
#   n <- length(y)
#   
#   ## Check dimension of P prior to calculating the basis
#   
#   if(n - dimBS(basis=basis,kernel=TRUE,degree=K[,1],segments=K[,2]) <= cv.df.min)
#     return(sqrt(.Machine$double.xmax))
#   
#   ## Otherwise, compute the cross-validation function
#   
#   if(is.null(z)) {
#     ## Here we need to use lm.wfit throughout with weights, but lm.wfit
#     ## needs non-null weights, so if weights are null create a vector of
#     ## ones
#     ## No categorical predictors, never reached when called by crs()
#     if(any(K[,1] > 0)) {
#       ## Check for rank-deficient fit
#       P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
#       ## Test for lack of degrees of freedom
#       if(NCOL(P) >= (n-1))
#         return(sqrt(.Machine$double.xmax))
#       if(basis=="additive" || basis=="glp") {
#         ## Test for full column rank
#         if(!singular.ok) {
#           if(!is.fullrank(cbind(1,P)))
#             return(sqrt(.Machine$double.xmax))
#         }
#         ## Additive spline regression models have an intercept in the lm()
#         ## model (though not in the gsl.bs function)
#         if(is.null(tau)) {
#           if(is.null(weights))
#             epsilon <- residuals(model <- lm.fit(cbind(1,P),y))
#           else
#             epsilon <- residuals(model <- lm.wfit(cbind(1,P),y,weights))
#         } else {
#           if(is.null(weights))
#             residuals <- tryCatch(residuals(rq.fit(cbind(1,P),y,tau=tau,method="fn")),error=function(e){FALSE})
#           else
#             residuals <- tryCatch(residuals(rq.wfit(cbind(1,P),y,tau=tau,weights,method="fn")),error=function(e){FALSE})
#           if(is.logical(residuals))
#             return(sqrt(.Machine$double.xmax))
#         }
#       } else {
#         ## Test for full column rank
#         if(!singular.ok) {
#           if(!is.fullrank(P))
#             return(sqrt(.Machine$double.xmax))
#         }
#         if(is.null(tau)) {
#           if(is.null(weights))
#             epsilon <- residuals(model <- lm.fit(P,y))
#           else
#             epsilon <- residuals(model <- lm.wfit(P,y,weights))
#         } else {
#           if(is.null(weights))
#             residuals <- tryCatch(residuals(rq.fit(P,y,tau=tau,method="fn")),error=function(e){FALSE})
#           else
#             residuals <- tryCatch(residuals(rq.wfit(P,y,tau=tau,weights,method="fn")),error=function(e){FALSE})
#           if(is.logical(residuals))
#             return(sqrt(.Machine$double.xmax))
#         }
#       }
#       htt <- hat(P)
#       htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
#     } else {
#       htt <- rep(1/n,n)
#       epsilon <- y-mean(y)
#     }
#     
#   } else {
#     
#     ## Categorical predictors - this is the workhorse
#     z <- as.matrix(z)
#     num.z <- NCOL(z)
#     epsilon <- numeric(length=n)
#     htt <- numeric(length=n)
#     ## At least one predictor for which degree > 0
#     if(any(K[,1] > 0)) {
#       P <- prod.spline(x=x,K=K,knots=knots,basis=basis,display.warnings=display.warnings)
#       ## Test for lack of degrees of freedom
#       if(NCOL(P) >= (n-1))
#         return(sqrt(.Machine$double.xmax))
#       
#       ## 2025: Hoist cbind out of the loop for efficiency
#       if(basis=="additive" || basis=="glp") {
#         XP <- cbind(1,P)
#       }
#       
#       for(i in 1:nrow.z.unique) {
#         zz <- ind == ind.vals[i]
#         L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
#         if(!is.null(weights)) L <- weights*L
#         if(basis=="additive" || basis=="glp") {
#           ## Test for full column rank
#           if(!singular.ok) {
#             ## Note: using XP (cbind(1,P))
#             if(!is.fullrank(XP*L))
#               return(sqrt(.Machine$double.xmax))
#           }
#           ## Additive spline regression models have an intercept in
#           ## the lm() model (though not in the gsl.bs function)
#           if(is.null(tau)) {
#             ## 2025: Optimize using .lm.fit with manual weighting
#             sw <- sqrt(L)
#             model <- .lm.fit(XP*sw,y*sw)
#             ## Check rank from model instead of is.fullrank
#             if(!singular.ok && model$rank < ncol(XP))
#               return(sqrt(.Machine$double.xmax))
#           } else {
#             ## Test for full column rank (rq case)
#             if(!singular.ok && !is.fullrank(XP*L))
#               return(sqrt(.Machine$double.xmax))
#             
#             model <- tryCatch(rq.wfit(XP,y,weights=L,tau=tau,method="fn"),error=function(e){FALSE})
#             if(is.logical(model))
#               return(sqrt(.Machine$double.xmax))
#             sw <- sqrt(L)
#             model.hat <- .lm.fit(XP*sw,y*sw)
#           }
#         } else {
#           if(is.null(tau)) {
#             sw <- sqrt(L)
#             model <- .lm.fit(P*sw,y*sw)
#             if(!singular.ok && model$rank < ncol(P))
#               return(sqrt(.Machine$double.xmax))
#           } else {
#             ## Test for full column rank (rq case)
#             if(!singular.ok && !is.fullrank(P*L))
#               return(sqrt(.Machine$double.xmax))
#             
#             model <- tryCatch(rq.wfit(P,y,weights=L,tau=tau,method="fn"),error=function(e){FALSE})
#             if(is.logical(model))
#               return(sqrt(.Machine$double.xmax))
#             sw <- sqrt(L)
#             model.hat <- .lm.fit(P*sw,y*sw)
#           }
#         }
#         
#         if(is.null(tau)) {
#           epsilon[zz] <- (model$residuals/sw)[zz]
#           htt[zz] <- hat.from.lm.fit(model)[zz]
#         } else {
#           epsilon[zz] <- residuals(model)[zz]
#           htt[zz] <- hat.from.lm.fit(model.hat)[zz]
#         }
#       }
#       
#     } else {
#       ## No predictors for which degree > 0
#       z.factor <- data.frame(factor(z[,1]),ordered=is.ordered.z[1])
#       if(num.z > 1) for(i in 2:num.z) z.factor <- data.frame(z.factor,factor(z[,i],ordered=is.ordered.z[i]))
#       
#       ## 2025: Hoist matrix creation out of loop
#       X0 <- matrix(1,n,1)
#       
#       for(i in 1:nrow.z.unique) {
#         zz <- ind == ind.vals[i]
#         L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
#         if(!is.null(weights)) L <- weights*L
#         
#         ## Whether we use additive, glp, or tensor products, this
#         ## model has no continuous predictors hence the intercept is
#         ## the parameter that may shift with the categorical
#         ## predictors
#         if(is.null(tau)) {
#           sw <- sqrt(L)
#           model <- .lm.fit(X0*sw,y*sw)
#           if(!singular.ok && model$rank < ncol(X0))
#             return(sqrt(.Machine$double.xmax))
#           htt[zz] <- hat.from.lm.fit(model)[zz]
#         } else {
#           ## Test for full column rank
#           if(!singular.ok && !is.fullrank(X0*L))
#             return(sqrt(.Machine$double.xmax))
#           
## The following function is a wrapper of cv.factor.spline to
## handle the situation when knots="auto". It will call
## cv.factor.spline two times and return the minimum with
## the optimal knots.

cv.factor.spline.wrapper <- function(x,
                                     y,
                                     z=NULL,
                                     K,
                                     I=NULL,
                                     knots=c("quantiles","uniform","auto"),
                                     basis=c("additive","tensor","glp"),
                                     cv.func=c("cv.ls","cv.gcv","cv.aic"),
                                     cv.df.min=1,
                                     tau=NULL,
                                     weights=NULL,
                                     singular.ok=FALSE,
                                     display.warnings=TRUE) {
  
  knots.opt <- knots
  
  if(knots == "auto") {
    
    knots.opt <- "quantiles"
    
    cv <- cv.factor.spline(x=x,
                           y=y,
                           z=z,
                           K=K,
                           I=I,
                           knots="quantiles",
                           basis=basis,
                           cv.func=cv.func,
                           cv.df.min=cv.df.min,
                           tau=tau,
                           weights=weights,
                           singular.ok=singular.ok,
                           display.warnings=display.warnings)
    
    cv.uniform <- cv.factor.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=K,
                                   I=I,
                                   knots="uniform",
                                   basis=basis,
                                   cv.func=cv.func,
                                   cv.df.min=cv.df.min,
                                   tau=tau,
                                   weights=weights,
                                   singular.ok=singular.ok,
                                   display.warnings=display.warnings)
    if(cv > cv.uniform) {
      cv <- cv.uniform
      knots.opt <- "uniform"
    }
    
  }	else {
    
    cv <- cv.factor.spline(x=x,
                           y=y,
                           z=z,
                           K=K,
                           I=I,
                           knots=knots,
                           basis=basis,
                           cv.func=cv.func,
                           cv.df.min=cv.df.min,
                           tau=tau,
                           weights=weights,
                           singular.ok=singular.ok,
                           display.warnings=display.warnings)
    
  }
  
  attr(cv, "knots.opt") <- knots.opt
  
  return(cv)
}

## Drop-in replacement for cv.factor.spline with improved handling of 
## rank-deficient and near-singular designs
##
## USAGE: Exactly the same as original cv.factor.spline, with additional
##        optional arguments after ... for controlling regularization
##
## NEW ARGUMENTS (all optional, appear after ...):
##   use.ridge           - Enable ridge regularization (default: TRUE)
##   ridge.lambda        - Ridge parameter (NULL = auto-select based on k/n ratio)
##   ridge.threshold     - Ratio k/n above which to use ridge (default: 0.7)
##   use.svd.fallback    - Use SVD if .lm.fit fails (default: TRUE)
##   smooth.penalty      - Use smooth penalties instead of hard cutoffs (default: TRUE)
##   penalty.scale       - Scale for smooth penalty (default: 1000)

cv.factor.spline <- function(x,
                             y,
                             z=NULL,
                             K,
                             I=NULL,
                             knots=c("quantiles","uniform"),
                             basis=c("additive","tensor","glp"),
                             cv.func=c("cv.ls","cv.gcv","cv.aic"),
                             cv.df.min=1,
                             tau=NULL,
                             weights=NULL,
                             singular.ok=FALSE,
                             display.warnings=TRUE,
                             ...,
                             use.ridge=TRUE,
                             ridge.lambda=NULL,
                             ridge.threshold=0.7,
                             use.svd.fallback=TRUE,
                             smooth.penalty=TRUE,
                             penalty.scale=1000) {
  
  if(missing(x) || missing(y) || missing(K)) stop(" must provide x, y and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  cv.func <- match.arg(cv.func)
  
  n <- NROW(x)
  have_tau <- !is.null(tau)
  have_w <- !is.null(weights)
  is_add <- (basis == "additive" || basis == "glp")
  cv.maxPenalty <- resolve_cv_maxPenalty(NULL, y, weights = weights)
  
  ## Check dimension of P prior to calculating the basis
  
  if(is.null(I)) {
    categories <- NULL
  } else {
    categories <- numeric()
    for(i in 1:NCOL(z)) categories[i] <- length(unique(z[,i]))
  }
  
  ## Calculate expected degrees of freedom
  k_expected <- dimBS(basis=basis, kernel=TRUE, degree=K[,1], 
                      segments=K[,2], include=I, categories=categories)
  
  ## IMPROVEMENT 1: Smooth penalty instead of hard cutoff
  df_remaining <- n - k_expected
  
  if(df_remaining <= cv.df.min) {
    if(smooth.penalty) {
      ## Smooth exponential penalty - allows optimizer to navigate
      deficit <- cv.df.min - df_remaining
      if(deficit > 10) {
        ## Completely infeasible
        return(cv.maxPenalty)
      } else {
        ## Smooth penalty that grows exponentially
        penalty_mult <- exp(deficit / 2)
        base_cv <- penalty.scale
        return(base_cv * penalty_mult)
      }
    } else {
      ## Original hard cutoff behavior
      return(cv.maxPenalty)
    }
  }
  
  ## Otherwise, compute the cross-validation function
  
  if(any(K[,1] > 0) || any(I > 0)) {
    P <- prod.spline(x=x, z=z, K=K, I=I, knots=knots, basis=basis, 
                     display.warnings=display.warnings)
    
    k_actual <- NCOL(P)
    
    ## IMPROVEMENT 2: More intelligent rank-deficiency check
    ## Check ratio k/n rather than just k >= n-1
    ratio <- k_actual / n
    
    if(ratio >= 0.99 && !use.ridge) {
      ## Extremely rank deficient, even ridge may not help
      if(smooth.penalty) {
        penalty_mult <- exp((ratio - 0.99) / 0.01)
        return(penalty.scale * penalty_mult)
      } else {
        return(cv.maxPenalty)
      }
    }
    
    ## Pre-calculate sw if weights exist
    if(have_w) sw <- sqrt(weights)
    
    ## Set up design matrix X
    if(is_add) {
      X <- cbind(1, P)
    } else {
      X <- P
    }
    
    ## IMPROVEMENT 3: Ridge regularization for near-singular cases
    use_ridge_now <- use.ridge && (ratio > ridge.threshold)
    
    if(use_ridge_now) {
      ## Auto-select lambda if not provided
      if(is.null(ridge.lambda)) {
        if(ratio > 0.95) {
          ridge.lambda <- 1e-2      ## Strong regularization
        } else if(ratio > 0.85) {
          ridge.lambda <- 1e-3      ## Moderate regularization
        } else if(ratio > 0.75) {
          ridge.lambda <- 5e-4      ## Light regularization
        } else {
          ridge.lambda <- 1e-4      ## Very light regularization
        }
      }
      
      ## Augment design matrix with ridge penalty
      ## Solve: (X'X + lambda*I)beta = X'y
      ## Equivalent to: [X; sqrt(lambda)I][beta] = [y; 0]
      k_X <- ncol(X)
      X_aug <- rbind(X, sqrt(ridge.lambda) * diag(k_X))
      y_aug <- c(y, rep(0, k_X))
      
      ## Apply weights to augmented system if needed
      if(have_w) {
        sw_aug <- c(sw, rep(1, k_X))
        X_fit <- X_aug * sw_aug
        y_fit <- y_aug * sw_aug
      } else {
        X_fit <- X_aug
        y_fit <- y_aug
      }
      
      n_aug <- length(y_fit)
      
    } else {
      ## No ridge regularization
      if(have_w) {
        X_fit <- X * sw
        y_fit <- y * sw
      } else {
        X_fit <- X
        y_fit <- y
      }
      n_aug <- n
    }
    
    ## IMPROVEMENT 4: Fit with error handling and SVD fallback
    if(!have_tau) {
      ## Least squares fitting
      
      model <- tryCatch({
        .lm.fit(X_fit, y_fit, tol=1e-7)
      }, error = function(e) {
        if(use.svd.fallback) {
          if(display.warnings) {
            warning("lm.fit failed, using SVD fallback")
          }
          ## SVD-based fitting
          svd_result <- svd_lm_fit(X_fit, y_fit, tol=1e-7)
          return(svd_result)
        } else {
          return(NULL)
        }
      })
      
      if(is.null(model)) {
        return(cv.maxPenalty)
      }
      
      ## Calculate residuals on ORIGINAL scale (not augmented)
      if(use_ridge_now) {
        ## Use original X, not augmented
        fitted_vals <- X %*% model$coefficients
        epsilon <- y - fitted_vals
      } else {
        if(have_w) {
          epsilon <- model$residuals / sw
        } else {
          epsilon <- model$residuals
        }
      }
      
      ## Check for rank deficiency (only if not using ridge)
      if(!singular.ok && !use_ridge_now) {
        if(model$rank < ncol(X)) {
          if(smooth.penalty) {
            ## Smooth penalty based on rank deficiency severity
            rank_deficit <- ncol(X) - model$rank
            penalty_mult <- exp(rank_deficit / ncol(X))
            return(penalty.scale * penalty_mult)
          } else {
            return(cv.maxPenalty)
          }
        }
      }
      
      ## Calculate hat values for CV
      ## For ridge regression, we need to adjust the hat matrix calculation
      if(use_ridge_now) {
        ## Effective hat matrix for ridge: H = X(X'X + lambda*I)^{-1}X'
        ## For augmented formulation: use hat values from augmented fit
        ## but only for original observations
        htt_aug <- hat.from.lm.fit(model)
        htt <- htt_aug[1:n]
      } else {
        htt <- hat.from.lm.fit(model)
      }
      
    } else {
      ## Quantile regression
      
      ## Check for rank deficiency
      if(!singular.ok && !is.fullrank(X))
        return(cv.maxPenalty)
      
      if(!have_w) {
        model <- tryCatch(
          rq.fit(X, y, tau=tau, method="fn"),
          error=function(e) {FALSE}
        )
        model.hat <- .lm.fit(X, y)
      } else {
        model <- tryCatch(
          rq.wfit(X, y, weights=weights, tau=tau, method="fn"),
          error=function(e) {FALSE}
        )
        model.hat <- .lm.fit(X*sw, y*sw)
      }
      
      if(is.logical(model))
        return(cv.maxPenalty)
      
      epsilon <- residuals(model)
      htt <- hat.from.lm.fit(model.hat)
    }
    
    ## Clamp hat values in-place
    idx <- htt >= 1
    if(any(idx)) htt[idx] <- 1 - .Machine$double.eps
    
  } else {
    ## No relevant predictors
    htt <- rep.int(1/n, n)
    epsilon <- y - mean(y)
  }
  
  ## If weights exist, need to use weighted residuals
  if(have_w) epsilon <- epsilon * sqrt(weights)
  
  ## Calculate cross-validation criterion
  if(cv.func == "cv.ls") {
    if(!have_tau)
      cv <- mean(epsilon^2 / (1 - htt)^2)
    else
      cv <- mean(check.function(epsilon, tau) / (1 - htt)^(1/sqrt(tau*(1-tau))))
  } else if(cv.func == "cv.gcv") {
    if(!have_tau)
      cv <- mean(epsilon^2 / (1 - mean(htt))^2)
    else
      cv <- mean(check.function(epsilon, tau) / (1 - mean(htt))^(1/sqrt(tau*(1-tau))))
  } else if(cv.func == "cv.aic") {
    traceH <- sum(htt)
    if(!have_tau) {
      sigmasq <- mean(epsilon^2)
      penalty <- ((1 + traceH/n) / (1 - (traceH + 2)/n))
    } else {
      sigmasq <- mean(check.function(epsilon, tau))
      penalty <- ((1 + traceH/n) / (1 - (traceH + 2)/n)) * (0.5/sqrt(tau*(1-tau)))
    }
    cv <- ifelse(penalty < 0, cv.maxPenalty, log(sigmasq) + penalty)
  }
  
  return(ifelse(!is.na(cv), cv, cv.maxPenalty))
}

## Drop-in replacement for cv.kernel.spline with improved handling of 
## rank-deficient and near-singular designs
##
## USAGE: Exactly the same as original cv.kernel.spline, with additional
##        optional arguments after ... for controlling regularization
##
## NEW ARGUMENTS (all optional, appear after ...):
##   use.ridge           - Enable ridge regularization (default: TRUE)
##   ridge.lambda        - Ridge parameter (NULL = auto-select based on k/n ratio)
##   ridge.threshold     - Ratio k/n above which to use ridge (default: 0.7)
##   use.svd.fallback    - Use SVD if .lm.fit fails (default: TRUE)
##   smooth.penalty      - Use smooth penalties instead of hard cutoffs (default: TRUE)
##   penalty.scale       - Scale for smooth penalty (default: 1000)

## Drop-in replacement for cv.kernel.spline with improved handling of 
## rank-deficient and near-singular designs
##
## USAGE: Exactly the same as original cv.kernel.spline, with additional
##        optional arguments after ... for controlling regularization
##
## NEW ARGUMENTS (all optional, appear after ...):
##   use.ridge           - Enable ridge regularization (default: TRUE)
##   ridge.lambda        - Ridge parameter (NULL = auto-select based on k/n ratio)
##   ridge.threshold     - Ratio k/n above which to use ridge (default: 0.7)
##   use.svd.fallback    - Use SVD if .lm.fit fails (default: TRUE)
##   smooth.penalty      - Use smooth penalties instead of hard cutoffs (default: TRUE)
##   penalty.scale       - Scale for smooth penalty (default: 1000)

cv.kernel.spline <- function(x,
                             y,
                             z=NULL,
                             K,
                             lambda=NULL,
                             z.unique,
                             ind,
                             ind.vals,
                             ind.list=NULL,
                             nrow.z.unique,
                             is.ordered.z=NULL,
                             knots=c("quantiles","uniform"),
                             basis=c("additive","tensor","glp"),
                             cv.func=c("cv.ls","cv.gcv","cv.aic"),
                             cv.df.min=1,
                             tau=NULL,
                             weights=NULL,
                             singular.ok=FALSE,
                             display.warnings=TRUE,
                             ...,
                             use.ridge=TRUE,
                             ridge.lambda=NULL,
                             ridge.threshold=0.7,
                             use.svd.fallback=TRUE,
                             smooth.penalty=TRUE,
                             penalty.scale=1000) {
  
  if(missing(x) || missing(y) || missing(K)) stop(" must provide x, y and K")
  
  if(!is.matrix(K)) stop(" K must be a two-column matrix")
  
  basis <- match.arg(basis)
  if(is.null(is.ordered.z)) stop(" is.ordered.z must be provided")
  knots <- match.arg(knots)
  cv.func <- match.arg(cv.func)
  
  ## Without computing P, compute the number of columns that P would
  ## be and if degrees of freedom is 1 or less, return a large penalty.
  
  n <- length(y)
  cv.maxPenalty <- resolve_cv_maxPenalty(NULL, y, weights = weights)
  
  ## Check dimension of P prior to calculating the basis
  k_expected <- dimBS(basis=basis, kernel=TRUE, degree=K[,1], segments=K[,2])
  
  ## IMPROVEMENT 1: Smooth penalty instead of hard cutoff
  df_remaining <- n - k_expected
  
  if(df_remaining <= cv.df.min) {
    if(smooth.penalty) {
      ## Smooth exponential penalty - allows optimizer to navigate
      deficit <- cv.df.min - df_remaining
      if(deficit > 10) {
        ## Completely infeasible
        return(cv.maxPenalty)
      } else {
        ## Smooth penalty that grows exponentially
        penalty_mult <- exp(deficit / 2)
        base_cv <- penalty.scale
        return(base_cv * penalty_mult)
      }
    } else {
      ## Original hard cutoff behavior
      return(cv.maxPenalty)
    }
  }
  
  ## Helper function for fitting with ridge/SVD
  fit_with_ridge <- function(X, y_fit, weights_fit, use_ridge_now, 
                             ridge_lambda, use_tau, sw=NULL) {
    
    k_X <- ncol(X)
    n_X <- nrow(X)
    
    if(use_ridge_now) {
      ## Augment design matrix
      X_aug <- rbind(X, sqrt(ridge_lambda) * diag(k_X))
      y_aug <- c(y_fit, rep(0, k_X))
      
      if(!is.null(sw)) {
        sw_aug <- c(sw, rep(1, k_X))
        X_fit <- X_aug * sw_aug
        y_fit_final <- y_aug * sw_aug
      } else {
        X_fit <- X_aug
        y_fit_final <- y_aug
      }
    } else {
      if(!is.null(sw)) {
        X_fit <- X * sw
        y_fit_final <- y_fit * sw
      } else {
        X_fit <- X
        y_fit_final <- y_fit
      }
    }
    
    ## Fit with error handling
    model <- tryCatch({
      .lm.fit(X_fit, y_fit_final, tol=1e-7)
    }, error = function(e) {
      if(use.svd.fallback) {
        if(display.warnings) {
          warning("lm.fit failed, using SVD fallback")
        }
        svd_result <- svd_lm_fit(X_fit, y_fit_final, tol=1e-7)
        return(svd_result)
      } else {
        return(NULL)
      }
    })
    
    return(model)
  }
  
  ## Otherwise, compute the cross-validation function
  
  if(is.null(z)) {
    ## No categorical predictors, never reached when called by crs()
    if(any(K[,1] > 0)) {
      P <- prod.spline(x=x, K=K, knots=knots, basis=basis, 
                       display.warnings=display.warnings)
      
      k_actual <- NCOL(P)
      ratio <- k_actual / n
      
      ## IMPROVEMENT 2: More intelligent rank-deficiency check
      if(ratio >= 0.99 && !use.ridge) {
        if(smooth.penalty) {
          penalty_mult <- exp((ratio - 0.99) / 0.01)
          return(penalty.scale * penalty_mult)
        } else {
          return(cv.maxPenalty)
        }
      }
      
      ## Determine if we should use ridge
      use_ridge_now <- use.ridge && (ratio > ridge.threshold)
      
      if(use_ridge_now) {
        ## Auto-select lambda if not provided
        if(is.null(ridge.lambda)) {
          if(ratio > 0.95) {
            ridge.lambda <- 1e-2
          } else if(ratio > 0.85) {
            ridge.lambda <- 1e-3
          } else if(ratio > 0.75) {
            ridge.lambda <- 5e-4
          } else {
            ridge.lambda <- 1e-4
          }
        }
      }
      
      if(basis=="additive" || basis=="glp") {
        X <- cbind(1, P)
        
        ## Test for full column rank (only if not using ridge)
        if(!singular.ok && !use_ridge_now) {
          if(!is.fullrank(X)) {
            if(smooth.penalty) {
              return(penalty.scale * 2)
            } else {
              return(cv.maxPenalty)
            }
          }
        }
        
        ## Additive spline regression models have an intercept
        if(is.null(tau)) {
          sw <- if(!is.null(weights)) sqrt(weights) else NULL
          model <- fit_with_ridge(X, y, weights, use_ridge_now, 
                                  ridge.lambda, FALSE, sw)
          
          if(is.null(model)) {
            return(cv.maxPenalty)
          }
          
          ## Calculate residuals on original scale
          if(use_ridge_now) {
            fitted_vals <- X %*% model$coefficients
            epsilon <- y - fitted_vals
          } else {
            if(!is.null(weights)) {
              epsilon <- model$residuals / sqrt(weights)
            } else {
              epsilon <- model$residuals
            }
          }
          
          ## Check rank (only if not using ridge)
          if(!singular.ok && !use_ridge_now) {
            if(model$rank < ncol(X)) {
              if(smooth.penalty) {
                rank_deficit <- ncol(X) - model$rank
                penalty_mult <- exp(rank_deficit / ncol(X))
                return(penalty.scale * penalty_mult)
              } else {
                return(cv.maxPenalty)
              }
            }
          }
          
        } else {
          ## Quantile regression case
          if(!is.null(weights))
            model <- tryCatch(rq.wfit(X, y, weights=weights, tau=tau, method="fn"),
                              error=function(e){FALSE})
          else
            model <- tryCatch(rq.fit(X, y, tau=tau, method="fn"),
                              error=function(e){FALSE})
          
          if(is.logical(model))
            return(cv.maxPenalty)
          
          epsilon <- residuals(model)
        }
        
      } else {
        ## Tensor basis
        X <- P
        
        if(!singular.ok && !use_ridge_now) {
          if(!is.fullrank(X)) {
            if(smooth.penalty) {
              return(penalty.scale * 2)
            } else {
              return(cv.maxPenalty)
            }
          }
        }
        
        if(is.null(tau)) {
          sw <- if(!is.null(weights)) sqrt(weights) else NULL
          model <- fit_with_ridge(X, y, weights, use_ridge_now, 
                                  ridge.lambda, FALSE, sw)
          
          if(is.null(model)) {
            return(cv.maxPenalty)
          }
          
          if(use_ridge_now) {
            fitted_vals <- X %*% model$coefficients
            epsilon <- y - fitted_vals
          } else {
            if(!is.null(weights)) {
              epsilon <- model$residuals / sqrt(weights)
            } else {
              epsilon <- model$residuals
            }
          }
          
          if(!singular.ok && !use_ridge_now) {
            if(model$rank < ncol(X)) {
              if(smooth.penalty) {
                rank_deficit <- ncol(X) - model$rank
                penalty_mult <- exp(rank_deficit / ncol(X))
                return(penalty.scale * penalty_mult)
              } else {
                return(cv.maxPenalty)
              }
            }
          }
          
        } else {
          if(!is.null(weights))
            model <- tryCatch(rq.wfit(X, y, weights=weights, tau=tau, method="fn"),
                              error=function(e){FALSE})
          else
            model <- tryCatch(rq.fit(X, y, tau=tau, method="fn"),
                              error=function(e){FALSE})
          
          if(is.logical(model))
            return(cv.maxPenalty)
          
          epsilon <- residuals(model)
        }
      }
      
      htt <- hat(P)
      htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
      
    } else {
      htt <- rep(1/n, n)
      epsilon <- y - mean(y)
    }
    
  } else {
    
    ## Categorical predictors - this is the workhorse
    if(!is.matrix(z)) z <- as.matrix(z)
    num.z <- NCOL(z)
    epsilon <- numeric(length=n)
    htt <- numeric(length=n)
    
    ## At least one predictor for which degree > 0
    if(any(K[,1] > 0)) {
      P <- prod.spline(x=x, K=K, knots=knots, basis=basis, 
                       display.warnings=display.warnings)
      
      k_actual <- NCOL(P)
      ratio <- k_actual / n
      
      ## IMPROVEMENT: Check ratio instead of hard threshold
      if(ratio >= 0.99 && !use.ridge) {
        if(smooth.penalty) {
          penalty_mult <- exp((ratio - 0.99) / 0.01)
          return(penalty.scale * penalty_mult)
        } else {
          return(cv.maxPenalty)
        }
      }
      
      ## Determine if we should use ridge
      use_ridge_now <- use.ridge && (ratio > ridge.threshold)
      
      if(use_ridge_now) {
        if(is.null(ridge.lambda)) {
          if(ratio > 0.95) {
            ridge.lambda <- 1e-2
          } else if(ratio > 0.85) {
            ridge.lambda <- 1e-3
          } else if(ratio > 0.75) {
            ridge.lambda <- 5e-4
          } else {
            ridge.lambda <- 1e-4
          }
        }
      }
      
      ## Hoist cbind out of the loop for efficiency
      if(basis=="additive" || basis=="glp") {
        XP <- cbind(1, P)
      } else {
        XP <- P
      }
      
      for(i in 1:nrow.z.unique) {
        if(!is.null(ind.list)) {
          zz <- ind.list[[i]]
        } else {
          zz <- ind == ind.vals[i]
        }
        L <- prod.kernel(Z=z, z=z.unique[ind.vals[i],], lambda=lambda, 
                         is.ordered.z=is.ordered.z)
        if(!is.null(weights)) L <- weights * L
        
        ## Calculate ratio for this subset
        ratio_subset <- ncol(XP) / sum(zz)
        use_ridge_subset <- use.ridge && (ratio_subset > ridge.threshold)
        
        if(use_ridge_subset && is.null(ridge.lambda)) {
          if(ratio_subset > 0.95) {
            ridge.lambda.subset <- 1e-2
          } else if(ratio_subset > 0.85) {
            ridge.lambda.subset <- 1e-3
          } else if(ratio_subset > 0.75) {
            ridge.lambda.subset <- 5e-4
          } else {
            ridge.lambda.subset <- 1e-4
          }
        } else {
          ridge.lambda.subset <- ridge.lambda
        }
        
        if(basis=="additive" || basis=="glp") {
          ## Test for full column rank (only if not using ridge)
          if(!singular.ok && !use_ridge_subset) {
            if(!is.fullrank(XP*L)) {
              if(smooth.penalty) {
                return(penalty.scale * 2)
              } else {
                return(cv.maxPenalty)
              }
            }
          }
          
          if(is.null(tau)) {
            sw <- sqrt(L)
            
            if(use_ridge_subset) {
              k_X <- ncol(XP)
              X_aug <- rbind(XP, sqrt(ridge.lambda.subset) * diag(k_X))
              y_aug <- c(y, rep(0, k_X))
              sw_aug <- c(sw, rep(1, k_X))
              
              model <- tryCatch({
                .lm.fit(X_aug*sw_aug, y_aug*sw_aug, tol=1e-7)
              }, error = function(e) {
                if(use.svd.fallback) {
                  svd_lm_fit(X_aug*sw_aug, y_aug*sw_aug, tol=1e-7)
                } else {
                  NULL
                }
              })
            } else {
              model <- tryCatch({
                .lm.fit(XP*sw, y*sw, tol=1e-7)
              }, error = function(e) {
                if(use.svd.fallback) {
                  svd_lm_fit(XP*sw, y*sw, tol=1e-7)
                } else {
                  NULL
                }
              })
            }
            
            if(is.null(model)) {
              return(cv.maxPenalty)
            }
            
            ## Check rank from model instead of is.fullrank
            if(!singular.ok && !use_ridge_subset && model$rank < ncol(XP)) {
              if(smooth.penalty) {
                rank_deficit <- ncol(XP) - model$rank
                penalty_mult <- exp(rank_deficit / ncol(XP))
                return(penalty.scale * penalty_mult)
              } else {
                return(cv.maxPenalty)
              }
            }
            
          } else {
            ## Test for full column rank (rq case)
            if(!singular.ok && !is.fullrank(XP*L)) {
              if(smooth.penalty) {
                return(penalty.scale * 2)
              } else {
                return(cv.maxPenalty)
              }
            }
            
            model <- tryCatch(rq.wfit(XP, y, weights=L, tau=tau, method="fn"),
                              error=function(e){FALSE})
            if(is.logical(model))
              return(cv.maxPenalty)
            sw <- sqrt(L)
            model.hat <- .lm.fit(XP*sw, y*sw)
          }
          
        } else {
          ## Tensor case
          if(is.null(tau)) {
            sw <- sqrt(L)
            
            if(use_ridge_subset) {
              k_X <- ncol(P)
              X_aug <- rbind(P, sqrt(ridge.lambda.subset) * diag(k_X))
              y_aug <- c(y, rep(0, k_X))
              sw_aug <- c(sw, rep(1, k_X))
              
              model <- tryCatch({
                .lm.fit(X_aug*sw_aug, y_aug*sw_aug, tol=1e-7)
              }, error = function(e) {
                if(use.svd.fallback) {
                  svd_lm_fit(X_aug*sw_aug, y_aug*sw_aug, tol=1e-7)
                } else {
                  NULL
                }
              })
            } else {
              model <- tryCatch({
                .lm.fit(P*sw, y*sw, tol=1e-7)
              }, error = function(e) {
                if(use.svd.fallback) {
                  svd_lm_fit(P*sw, y*sw, tol=1e-7)
                } else {
                  NULL
                }
              })
            }
            
            if(is.null(model)) {
              return(cv.maxPenalty)
            }
            
            if(!singular.ok && !use_ridge_subset && model$rank < ncol(P)) {
              if(smooth.penalty) {
                rank_deficit <- ncol(P) - model$rank
                penalty_mult <- exp(rank_deficit / ncol(P))
                return(penalty.scale * penalty_mult)
              } else {
                return(cv.maxPenalty)
              }
            }
            
          } else {
            ## Test for full column rank (rq case)
            if(!singular.ok && !is.fullrank(P*L)) {
              if(smooth.penalty) {
                return(penalty.scale * 2)
              } else {
                return(cv.maxPenalty)
              }
            }
            
            model <- tryCatch(rq.wfit(P, y, weights=L, tau=tau, method="fn"),
                              error=function(e){FALSE})
            if(is.logical(model))
              return(cv.maxPenalty)
            sw <- sqrt(L)
            model.hat <- .lm.fit(P*sw, y*sw)
          }
        }
        
        if(is.null(tau)) {
          if(use_ridge_subset) {
            ## For ridge, calculate residuals on original scale
            if(basis=="additive" || basis=="glp") {
              fitted_vals <- XP %*% model$coefficients
            } else {
              fitted_vals <- P %*% model$coefficients
            }
            epsilon[zz] <- (y - fitted_vals)[zz]
            ## For ridge: hat values from augmented system, take first n
            htt_all <- hat.from.lm.fit(model)
            htt[zz] <- htt_all[1:n][zz]
          } else {
            epsilon[zz] <- (model$residuals/sw)[zz]
            htt[zz] <- hat.from.lm.fit(model)[zz]
          }
        } else {
          epsilon[zz] <- residuals(model)[zz]
          htt[zz] <- hat.from.lm.fit(model.hat)[zz]
        }
      }
      
    } else {
      ## No predictors for which degree > 0
      z.factor <- data.frame(factor(z[,1]), ordered=is.ordered.z[1])
      if(num.z > 1) for(i in 2:num.z) 
        z.factor <- data.frame(z.factor, factor(z[,i], ordered=is.ordered.z[i]))
      
      ## Hoist matrix creation out of loop
      X0 <- matrix(1, n, 1)
      
      for(i in 1:nrow.z.unique) {
        if(!is.null(ind.list)) {
          zz <- ind.list[[i]]
        } else {
          zz <- ind == ind.vals[i]
        }
        L <- prod.kernel(Z=z, z=z.unique[ind.vals[i],], lambda=lambda, 
                         is.ordered.z=is.ordered.z)
        if(!is.null(weights)) L <- weights * L
        
        if(is.null(tau)) {
          sw <- sqrt(L)
          model <- .lm.fit(X0*sw, y*sw)
          if(!singular.ok && model$rank < ncol(X0))
            return(cv.maxPenalty)
          htt[zz] <- hat.from.lm.fit(model)[zz]
        } else {
          ## Test for full column rank
          if(!singular.ok && !is.fullrank(X0*L))
            return(cv.maxPenalty)
          
          model <- tryCatch(rq.wfit(X0, y, weights=L, tau=tau, method="fn"),
                            error=function(e){FALSE})
          if(is.logical(model))
            return(cv.maxPenalty)
          sw <- sqrt(L)
          model.hat <- .lm.fit(X0*sw, y*sw)
          htt[zz] <- hat.from.lm.fit(model.hat)[zz]
        }
        if(is.null(tau))
          epsilon[zz] <- (model$residuals/sw)[zz]
        else
          epsilon[zz] <- residuals(model)[zz]
      }
    }
    
    htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
    
  }
  
  ## If weights exist, need to use weighted residuals
  if(!is.null(weights)) epsilon <- epsilon * sqrt(weights)
  
  if(cv.func == "cv.ls") {
    if(is.null(tau))
      cv <- mean(epsilon^2 / (1 - htt)^2)
    else
      cv <- mean(check.function(epsilon, tau) / (1 - htt)^(1/sqrt(tau*(1-tau))))
  } else if(cv.func == "cv.gcv") {
    if(is.null(tau))
      cv <- mean(epsilon^2 / (1 - mean(htt))^2)
    else
      cv <- mean(check.function(epsilon, tau) / (1 - mean(htt))^(1/sqrt(tau*(1-tau))))
  } else if(cv.func == "cv.aic") {
    traceH <- sum(htt)
    if(is.null(tau)) {
      sigmasq <- mean(epsilon^2)
      penalty <- ((1 + traceH/n) / (1 - (traceH + 2)/n))
    } else {
      sigmasq <- mean(check.function(epsilon, tau))
      penalty <- ((1 + traceH/n) / (1 - (traceH + 2)/n)) * (0.5/sqrt(tau*(1-tau)))
    }
    cv <- ifelse(penalty < 0, cv.maxPenalty, log(sigmasq) + penalty)
  }
  
  return(ifelse(!is.na(cv), cv, cv.maxPenalty))
}

## ============================================================================
## Helper function: SVD-based least squares (for fallback)
## ============================================================================

svd_lm_fit <- function(x, y, tol = 1e-7) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  ## Compute SVD: X = UDV'
  svd_x <- svd(x)
  
  ## Determine rank based on tolerance
  d <- svd_x$d
  rank <- sum(d > max(tol * d[1], 0))
  
  if(rank == 0) {
    ## Completely rank deficient
    return(list(
      coefficients = rep(0, p),
      residuals = y,
      fitted.values = rep(0, n),
      rank = 0,
      qr = matrix(0, n, p),
      qraux = rep(0, p),
      pivot = 1:p,
      tol = tol,
      effects = rep(0, n)
    ))
  }
  
  ## Truncate to effective rank
  d_inv <- ifelse(d > max(tol * d[1], 0), 1/d, 0)
  
  ## Compute coefficients: beta = V D^{-1} U' y
  u_truncated <- svd_x$u[, 1:rank, drop = FALSE]
  v_truncated <- svd_x$v[, 1:rank, drop = FALSE]
  
  coefficients <- v_truncated %*% (d_inv[1:rank] * (t(u_truncated) %*% y))
  
  ## Ensure full length coefficient vector
  if(length(coefficients) < p) {
    coef_full <- rep(0, p)
    coef_full[1:length(coefficients)] <- coefficients
    coefficients <- coef_full
  }
  
  fitted.values <- x %*% coefficients
  residuals <- y - fitted.values
  
  ## Create QR-like output for compatibility with hat.from.lm.fit
  ## We need to provide qr, qraux, pivot, tol, rank
  ## For SVD, we can construct an equivalent QR representation
  
  ## Create a mock QR object that will work with hat.from.lm.fit
  ## The key is that qr() should give us the right hat matrix
  qr_obj <- qr(x, tol = tol)
  
  return(list(
    coefficients = as.vector(coefficients),
    residuals = as.vector(residuals),
    fitted.values = as.vector(fitted.values),
    rank = rank,
    qr = qr_obj$qr,
    qraux = qr_obj$qraux,
    pivot = qr_obj$pivot,
    tol = tol,
    effects = d[1:min(rank, length(d))]
  ))
}
