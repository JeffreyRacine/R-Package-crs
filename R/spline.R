## We support multivariate numeric and categorical
## predictors. Multivariate splines are additive and tensor. Can be
## used for training or evaluation data. For continuous datatypes uses
## the B-spline, for factor datatypes uses indicator splines obtained
## from model.matrix(). It also computes derivatives for the
## continuous variables of arbitrary order (issues warning when order
## exceeds degree of spline) with interaction if specified.

## Complexity can be modified via the number of knots (segments) or the
## spline degree (degree)

prod.spline <- function(x,
                        z=NULL,
                        K=NULL,
                        I=NULL,
                        xeval=NULL,
                        zeval=NULL,
                        knots=c("quantiles","uniform"),
                        basis=c("additive","tensor","glp"),
                        deriv.index=1,
                        deriv=0) {

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
  if(deriv > K[deriv.index,1]) warning(" deriv order too large, result will be zero")
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

predict.kernel.spline <- function(x,
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
                                  tau=NULL){

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
  console <- printPush("Working...",console = console)

  model <- NULL ## Returned if model=FALSE and there exist categorical
                ## predictors

  if(is.null(z)) {

    ## First no categorical predictor case

    if(any(K[,1] > 0)) {

      ## Degree > 0

      P <- prod.spline(x=x,K=K,knots=knots,basis=basis)

      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))
          model <- lm(y~P)
        else
          model <- rq(y~P,tau=tau)          
      } else {
        if(is.null(tau))        
          model <- lm(y~P-1)
        else
          model <- rq(y~P-1,tau=tau)          
      }
      if(is.null(xeval)) {
        fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
      } else {
        P <- prod.spline(x=x,K=K,xeval=xeval,knots=knots,basis=basis)
        fit.spline <- predict(model,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
      }

    } else {

      ## Degree == 0

      if(is.null(tau))        
        model <- lm(y~1)
      else 
        model <- rq(y~1,tau=tau)      
      if(is.null(xeval)) {
        fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
      } else {
        fit.spline <- predict(model,newdata=data.frame(rep(coef(model),NROW(xeval))),interval="confidence",se.fit=TRUE)
      }
    }

    fit.spline <- cbind(fit.spline[[1]],se=fit.spline[[2]])

    if(is.null(tau))        
      htt <- hatvalues(model)
    else
      htt <- hat(model$x)
    
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
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis)
          k <- NCOL(P)
          if(basis=="additive" || basis=="glp") {
            if(is.null(tau))        
              model.z.unique <- lm(y~P,weights=L)
            else
              model.z.unique <- rq(y~P,weights=L,tau=tau)              
              model.z.unique.hat <- lm(y~P,weights=L)              
          } else {
            if(is.null(tau))        
              model.z.unique <- lm(y~P-1,weights=L)
            else
              model.z.unique <- rq(y~P-1,weights=L,tau=tau)              
              model.z.unique.hat <- lm(y~P-1,weights=L)              
          }
          if(model.return) model[[i]] <- model.z.unique
          if(is.null(tau))        
            htt[zz] <- hatvalues(model.z.unique)[zz]
          else
            htt[zz] <- hatvalues(model.z.unique.hat)[zz]
          
          P.hat[zz] <- sum(L)
          P <- prod.spline(x=x,K=K,xeval=x[zz,,drop=FALSE],knots=knots,basis=basis)
          tmp <- predict(model.z.unique,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
          fit.spline[zz,] <- cbind(tmp[[1]],tmp[[2]])
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
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis)
          k <- NCOL(P)
          if(basis=="additive" || basis=="glp") {
            if(is.null(tau))        
              model.z.unique <- lm(y~P,weights=L)
            else
              model.z.unique <- rq(y~P,weights=L,tau=tau)              
          } else {
            if(is.null(tau))        
              model.z.unique <- lm(y~P-1,weights=L)
            else
              model.z.unique <- rq(y~P-1,weights=L,tau=tau)              
          }
          if(model.return) model[[i]] <- model.z.unique
          P <- prod.spline(x=x,K=K,xeval=xeval[zz,,drop=FALSE],knots=knots,basis=basis)
          tmp <- predict(model.z.unique,newdata=data.frame(as.matrix(P)),interval="confidence",se.fit=TRUE)
          fit.spline[zz,] <- cbind(tmp[[1]],tmp[[2]])
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
          k <- 0
          ## Whether we use additive, glp, or tensor products, this
          ## model has no continuous predictors hence the intercept is
          ## the parameter that may shift with the categorical
          ## predictors
          if(is.null(tau))        
            model.z.unique <- lm(y~x.intercept-1,weights=L)
          else
            model.z.unique <- rq(y~x.intercept-1,weights=L,tau=tau)
            model.z.unique.hat <- lm(y~x.intercept-1,weights=L)
          if(model.return) model[[i]] <- model.z.unique
          if(is.null(tau))        
            htt[zz] <- hatvalues(model.z.unique)[zz]
          else
            htt[zz] <- hatvalues(model.z.unique.hat)[zz]
          P.hat[zz] <- sum(L)
          tmp <- predict(model.z.unique,newdata=data.frame(x.intercept=x.intercept[zz]),interval="confidence",se.fit=TRUE)
          fit.spline[zz,] <- cbind(tmp[[1]],tmp[[2]])
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
          k <- 0
          if(is.null(tau))        
            model.z.unique <- lm(y~x.intercept-1,weights=L)
          else
            model.z.unique <- rq(y~x.intercept-1,weights=L,tau=tau)            
          if(model.return) model[[i]] <- model.z.unique
          tmp <- predict(model.z.unique,newdata=data.frame(x.intercept=rep(1,num.eval)[zz]),interval="confidence",se.fit=TRUE)
          fit.spline[zz,] <- cbind(tmp[[1]],as.matrix(tmp[[2]]))
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

deriv.kernel.spline <- function(x,
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
                                tau=NULL) {

  if(deriv == 0) stop(" deriv must be greater than zero")

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")  

  basis <- match.arg(basis)
  knots <- match.arg(knots)
  if(is.null(is.ordered.z)) stop(" is.ordered.z must be provided")

  x <- as.matrix(x)

  ## Univariate additive spline bases have one less column than
  ## univariate tensor spline bases. This is used only for setting
  ## appropriate columns for derivative computation.
  
  if(basis=="additive" || basis=="glp") {
    K.additive <- K
    K.additive[,1] <- ifelse(K[,1]>0,K[,1]-1,K[,1])
  }

  if(!is.null(z)) z <- as.matrix(z)

  if(is.null(z)) {

    ## First no categorical predictor case

    if(K[deriv.index,1]!=0) {

      P <- prod.spline(x=x,K=K,knots=knots,basis=basis)      
      P.deriv <- prod.spline(x=x,K=K,xeval=xeval,knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv)
      dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
      dim.P.tensor <- NCOL(P)

      if(basis=="additive") {
        if(is.null(tau))
          model <- lm(y~P)
        else
          model <- rq(y~P,tau=tau)
          
        dim.P.deriv <- sum(K.additive[deriv.index,])
        deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)      
        deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
        deriv.ind.vec <- deriv.start:deriv.end
        deriv.spline <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]

        if(is.null(tau))        
          vcov.model <- vcov(model)[-1,-1,drop=FALSE]
        else
          vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]
        
        se.deriv <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
      } else if(basis=="tensor") {
        if(is.null(tau))
          model <- lm(y~P-1)
        else
          model <- rq(y~P-1,tau=tau)
          
        deriv.spline <- P.deriv%*%coef(model)

        if(is.null(tau))        
          vcov.model <- vcov(model)
        else
          vcov.model <- summary(model,covariance=TRUE)$cov

        se.deriv <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model.%*%t(P.deriv[i,,drop=FALSE])) })
      } else if(basis=="glp") {
        if(is.null(tau))
          model <- lm(y~P)
        else
          model <- rq(y~P,tau=tau)          
        deriv.spline <- P.deriv%*%coef(model)[-1]
        
        if(is.null(tau))        
          vcov.model <- vcov(model)[-1,-1,drop=FALSE]
        else
          vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]

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
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis)          
          P.deriv <- prod.spline(x=x,K=K,xeval=x[zz,,drop=FALSE],knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv)
          k <- NCOL(P)
          dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
          dim.P.tensor <- NCOL(P)

          if(basis=="additive") {
            if(is.null(tau))
              model <- lm(y~P)
            else
              model <- rq(y~P,tau=tau)              
            dim.P.deriv <- sum(K.additive[deriv.index,])
            deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)      
            deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
            deriv.ind.vec <- deriv.start:deriv.end
            deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]

            if(is.null(tau))        
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]

            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })            
          } else if(basis=="tensor") {
            if(is.null(tau))
              model <- lm(y~P-1)
            else
              model <- rq(y~P-1,tau=tau)            

            deriv.spline[zz] <- P.deriv%*%coef(model)

            if(is.null(tau))        
              vcov.model <- vcov(model)
            else
              vcov.model <- summary(model,covariance=TRUE)$cov
            
            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
          } else if(basis=="glp") {
            if(is.null(tau))
              model <- lm(y~P)
            else
              model <- rq(y~P,tau=tau)
              
            deriv.spline[zz] <- P.deriv%*%coef(model)[-1]

            if(is.null(tau))        
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]

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
          P <- prod.spline(x=x,K=K,knots=knots,basis=basis)
          P.deriv <- prod.spline(x=x,K=K,xeval=xeval[zz,,drop=FALSE],knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv)
          k <- NCOL(P)
          dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
          dim.P.tensor <- NCOL(P)

          if(basis=="additive") {
            if(is.null(tau))
              model <- lm(y~P,weights=L)
            else
              model <- rq(y~P,weights=L,tau=tau)              

            dim.P.deriv <- sum(K.additive[deriv.index,])
            deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)      
            deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
            deriv.ind.vec <- deriv.start:deriv.end
            deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]

            if(is.null(tau))        
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]

            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,deriv.ind.vec,drop=FALSE]%*%vcov.model[deriv.ind.vec,deriv.ind.vec]%*%t(P.deriv[i,deriv.ind.vec,drop=FALSE])) })
          } else if(basis=="tensor") {
            if(is.null(tau))
              model <- lm(y~P-1,weights=L)
            else
              model <- rq(y~P-1,weights=L,tau=tau)              
            deriv.spline[zz] <- P.deriv%*%coef(model)

            if(is.null(tau))        
              vcov.model <- vcov(model)
            else
              vcov.model <- summary(model,covariance=TRUE)$cov

            se.deriv[zz] <- sapply(1:NROW(P.deriv), function(i){ sqrt(P.deriv[i,,drop=FALSE]%*%vcov.model%*%t(P.deriv[i,,drop=FALSE])) })
          } else if(basis=="glp") {
            if(is.null(tau))
              model <- lm(y~P,weights=L)
            else
              model <- rq(y~P,weights=L,tau=tau)
              
            deriv.spline[zz] <- P.deriv%*%coef(model)[-1]

            if(is.null(tau))        
              vcov.model <- vcov(model)[-1,-1,drop=FALSE]
            else
              vcov.model <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]

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

predict.factor.spline <- function(x,
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
                                  tau=NULL){

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
  console <- printPush("Working...",console = console)

  if(any(K[,1] > 0)||any(I>0)) {

    ## Degree > 0

    P <- prod.spline(x=x,z=z,K=K,I=I,knots=knots,basis=basis)

    if(prune && is.null(prune.index)) {

      ## Pruning via step-wise CV but returning the pruned model only
      ## if the cross-validation score is improved (lower). We create
      ## a data frame so that we can readily determine columns that
      ## have been removed and assign logical values to all columns in
      ## P.
      P.df <- data.frame(P)
      names(P.df) <- paste("P",seq(1,NCOL(P.df)),sep="")
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau)) 
          model <- lm(y~.,data=P.df)
        else
          model <- rq(y~.,data=P.df,tau=tau)          
      } else {
        if(is.null(tau)) 
          model <- lm(y~.-1,data=P.df)
        else
          model <- rq(y~.-1,data=P.df,tau=tau)          
      }
      if(is.null(tau)) 
        cv <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      else
        cv <- cv.rq(model,tau=tau)
      console <- printClear(console)
      console <- printPush("Pruning...",console = console)
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))
          model.pruned <- stepCV(lm(y~.,data=P.df),
                                 scope=list(upper=~.,lower=~1),
                                 k=log(length(y)),
                                 trace=trace)
        else
          model.pruned <- stepCV(rq(y~.,data=P.df,tau=tau),
                                 scope=list(upper=~.,lower=~1),
                                 k=log(length(y)),
                                 trace=trace)

          
      } else {
        if(is.null(tau))
          model.pruned <- stepCV(lm(y~.-1,data=P.df),
                                 scope=list(upper=~.,lower=~1),
                                 k=log(length(y)),
                                 trace=trace)
        else
          model.pruned <- stepCV(rq(y~.-1,data=P.df,tau=tau),
                                 scope=list(upper=~.,lower=~1),
                                 k=log(length(y)),
                                 trace=trace)
          
      }
      if(is.null(tau))
        cv.pruned <- mean(residuals(model.pruned)^2/(1-hatvalues(model.pruned))^2)
      else
        cv.pruned <- cv.rq(model.pruned,tau=tau)
      
      if(cv.pruned <= cv) {
        IND <- logical()
        for(i in 1:NCOL(P.df)) IND[i] <- any(names(P.df)[i]==names(model.pruned$model[,-1,drop=FALSE]))
        if(basis=="additive" || basis=="glp") {
          if(is.null(tau))
            model <- lm(y~P[,IND,drop=FALSE])
          else
            model <- rq(y~P[,IND,drop=FALSE],tau=tau)            
        } else {
          if(is.null(tau))          
            model <- lm(y~P[,IND,drop=FALSE]-1)
          else
            model <- rq(y~P[,IND,drop=FALSE]-1,tau=tau)            
        }
      } else {
        warning(" pruned model did not lower cross-validation score, using non-pruned bases")
        IND <- !logical(length=NCOL(P))
        if(basis=="additive" || basis=="glp") {
          if(is.null(tau))          
            model <- lm(y~P)
          else
            model <- rq(y~P,tau=tau)            
        } else {
          if(is.null(tau))          
            model <- lm(y~P-1)
          else
            model <- rq(y~P-1,tau=tau)            
        }
      }
    } else if(prune) {
      ## Pruning, index passed in...
      IND <- prune.index
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))          
          model <- lm(y~P[,IND,drop=FALSE])
        else
          model <- rq(y~P[,IND,drop=FALSE],tau=tau)          
      } else {
        if(is.null(tau))          
          model <- lm(y~P[,IND,drop=FALSE]-1)
        else
          model <- rq(y~P[,IND,drop=FALSE]-1,tau=tau)          
      }
      cv <- NULL
      if(is.null(tau))          
        cv.pruned <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      else
        cv.pruned <- cv.rq(model,tau=tau)
    } else {
      ## No pruning
      IND <- !logical(length=NCOL(P))
      if(basis=="additive" || basis=="glp") {
        if(is.null(tau))          
          model <- lm(y~P)
        else
          model <- rq(y~P,tau=tau)          
      } else {
        if(is.null(tau))          
          model <- lm(y~P-1)
        else
          model <- rq(y~P-1,tau=tau)          
      }
      cv.pruned <- NULL
      if(is.null(tau))          
        cv <- mean(residuals(model)^2/(1-hatvalues(model))^2)
      else
        cv <- cv.rq(model,tau=tau)
    }      

    if(is.null(xeval)) {
      fit.spline <- predict(model,interval="confidence",se.fit=TRUE)
    } else {
      P <- prod.spline(x=x,z=z,K=K,I=I,xeval=xeval,zeval=zeval,knots=knots,basis=basis)
      fit.spline <- predict(model,newdata=data.frame(as.matrix(P[,IND,drop=FALSE])),interval="confidence",se.fit=TRUE)
    }

  } else {

    ## Degree == 0, no pruning possible
    IND <- TRUE

    if(is.null(tau))          
      model <- lm(y~1)
    else
      model <- rq(y~1,tau=tau)

    
    if(is.null(tau))
      cv <- mean(residuals(model)^2/(1-hatvalues(model))^2) ## Added
    else
      cv <- cv.rq(model,tau=tau)
      
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
    fit.spline <- cbind(fit.spline,NA)

  console <- printClear(console)
  console <- printPop(console)

  if(is.null(tau)) 
    htt <- hatvalues(model)
  else
    htt <- hat(model$x)

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

deriv.factor.spline <- function(x,
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
                                tau=NULL) {

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(deriv == 0) stop(" derivative must be a positive integer")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")  

  basis <- match.arg(basis)
  knots <- match.arg(knots)

  x <- as.matrix(x)

  ## Univariate additive spline bases have one less column than
  ## univariate tensor spline bases. This is used only for setting
  ## appropriate columns for derivative computation.
  
  if(basis=="additive" || basis=="glp") {
    K.additive <- K
    K.additive[,1] <- ifelse(K[,1]>0,K[,1]-1,K[,1])
  }

  if(K[deriv.index,1]!=0) {

    ## Degree > 0

    ## Estimate model on training data.
    
    P <- prod.spline(x=x,z=z,K=K,I=I,knots=knots,basis=basis)    
    P.deriv <- prod.spline(x=x,z=z,K=K,I=I,xeval=xeval,zeval=zeval,knots=knots,basis=basis,deriv.index=deriv.index,deriv=deriv)

    dim.P.no.tensor <- attr(P.deriv,"dim.P.no.tensor")
    dim.P.tensor <- NCOL(P)
    deriv.ind.vec <- logical(length=NCOL(P)) ## All false

    if(is.null(prune.index)) prune.index <- !logical(NCOL(P))

    ## Pad the following for proper handling of pruning

    coef.vec.model <- numeric(length=NCOL(P))
    vcov.mat.model <- matrix(0,nrow=NCOL(P),ncol=NCOL(P))

    if(basis=="additive") {
      if(is.null(tau))
        model <- lm(y~P[,prune.index,drop=FALSE])
      else
        model <- rq(y~P[,prune.index,drop=FALSE],tau=tau)        

      coef.vec.model[prune.index] <- coef(model)[-1]

      if(is.null(tau))        
        vcov.mat.model[prune.index,prune.index] <- vcov(model)[-1,-1,drop=FALSE]
      else
        vcov.mat.model[prune.index,prune.index] <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]

      dim.P.deriv <- sum(K.additive[deriv.index,])
      deriv.start <- ifelse(deriv.index!=1,sum(K.additive[1:(deriv.index-1),])+1,1)      
      deriv.end <- deriv.start+sum(K.additive[deriv.index,])-1
      deriv.ind.vec[deriv.start:deriv.end] <- TRUE
      deriv.ind.vec <- ifelse(prune.index,deriv.ind.vec,FALSE)
    } else if(basis=="tensor") {
      if(is.null(tau))
        model <- lm(y~P[,prune.index,drop=FALSE]-1)
      else
        model <- rq(y~P[,prune.index,drop=FALSE]-1,tau=tau)        
      coef.vec.model[prune.index] <- coef(model)
      
      if(is.null(tau))        
        vcov.mat.model[prune.index,prune.index] <- vcov(model)
      else
        vcov.mat.model[prune.index,prune.index] <- summary(model,covariance=TRUE)$cov

      deriv.ind.vec[1:dim.P.tensor] <- TRUE            
      deriv.ind.vec <- ifelse(prune.index,deriv.ind.vec,FALSE)
    } else if(basis=="glp") {
      if(is.null(tau))
        model <- lm(y~P[,prune.index,drop=FALSE])
      else
        model <- rq(y~P[,prune.index,drop=FALSE],tau=tau)        
      coef.vec.model[prune.index] <- coef(model)[-1]

      if(is.null(tau))        
        vcov.mat.model[prune.index,prune.index] <- vcov(model)[-1,-1,drop=FALSE]
      else
        vcov.mat.model[prune.index,prune.index] <- summary(model,covariance=TRUE)$cov[-1,-1,drop=FALSE]

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
                                     nrow.z.unique,
                                     is.ordered.z=NULL,
                                     knots=c("quantiles","uniform","auto"),
                                     basis=c("additive","tensor","glp"),
                                     cv.func=c("cv.ls","cv.gcv","cv.aic","cv.rq"),
                                     tau=NULL) {

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
                           nrow.z.unique=nrow.z.unique, 
                           is.ordered.z=is.ordered.z, 
                           knots="quantiles", 
                           basis=basis, 
                           cv.func=cv.func,
                           tau=tau)
    
    cv.uniform <- cv.kernel.spline(x=x, 
                                   y=y, 
                                   z=z, 
                                   K=K, 
                                   lambda=lambda, 
                                   z.unique=z.unique, 
                                   ind=ind, 
                                   ind.vals=ind.vals, 
                                   nrow.z.unique=nrow.z.unique, 
                                   is.ordered.z=is.ordered.z, 
                                   knots="uniform", 
                                   basis=basis, 
                                   cv.func=cv.func,
                                   tau=tau)
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
                           nrow.z.unique=nrow.z.unique, 
                           is.ordered.z=is.ordered.z, 
                           knots=knots, 
                           basis=basis, 
                           cv.func=cv.func,
                           tau=tau)

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

cv.kernel.spline <- function(x,
														 y,
														 z=NULL,
														 K,
														 lambda=NULL,
														 z.unique,
														 ind,
														 ind.vals,
														 nrow.z.unique,
														 is.ordered.z=NULL,
														 knots=c("quantiles","uniform"),
														 basis=c("additive","tensor","glp"),
														 cv.func=c("cv.ls","cv.gcv","cv.aic","cv.rq"),
                             tau=NULL) {

  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  
  if(!is.matrix(K)) stop(" K must be a two-column matrix")  
  
  basis <- match.arg(basis)
  if(is.null(is.ordered.z)) stop(" is.ordered.z must be provided")
  knots <- match.arg(knots)
  cv.func <- match.arg(cv.func)
  
  ## Without computing P, compute the number of columns that P would
  ## be and if degrees of freedom is 1 or less, return a large penalty.
  
  n <- length(y)
  
  ## Otherwise, compute the cross-validation function
  
  if(is.null(z)) {
    ## No categorical predictors
    if(any(K[,1] > 0)) {
      P <- prod.spline(x=x,K=K,knots=knots,basis=basis)
      if(NCOL(P) >= (n-1)) return(sqrt(.Machine$double.xmax))
      
      if(basis=="additive" || basis=="glp") {
        ## Additive spline regression models have an intercept in the lm()
        ## model (though not in the gsl.bs function)
        if(is.null(tau)) {
          epsilon <- residuals(lm.fit(cbind(1,P),y))
          ## Check for rank-deficient fit
          if(model$rank != (NCOL(P)+1)) return(sqrt(.Machine$double.xmax))
        } else {
          epsilon <- residuals(rq.fit(cbind(1,P),y,tau=tau))
          ## Check for rank-deficient fit                
          if(NCOL(model$x) != (NCOL(P)+1)) return(sqrt(.Machine$double.xmax))
        }
      } else {
        if(is.null(tau)) {
          epsilon <- residuals(lm.fit(P,y))
          ## Check for rank-deficient fit
          if(model$rank != NCOL(P)) return(sqrt(.Machine$double.xmax))
        } else {
          epsilon <- residuals(rq.fit(P,y,tau=tau))
          ## Check for rank-deficient fit                
          if(NCOL(model$x) != NCOL(P)) return(sqrt(.Machine$double.xmax))
        }
      }
      htt <- hat(P)
      htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
    } else {
      htt <- rep(1/n,n)
      epsilon <- y-mean(y)
    }
    if(cv.func == "cv.ls") {
      cv <- mean(epsilon^2/(1-htt)^2)
    } else if(cv.func == "cv.gcv"){
      cv <- mean(epsilon^2/(1-mean(htt))^2)
    } else if(cv.func == "cv.aic"){
      sigmasq <- mean(epsilon^2)
      traceH <- sum(htt)
      penalty <- (1+traceH/n)/(1-(traceH+2)/n)
      cv <- ifelse(penalty < 0, .Machine$double.xmax, log(sigmasq)+penalty);
    } else if(cv.func == "cv.rq") {
      ## Note - this is defined in util.R so if you modify there you must modify here also
      cv <- mean(check.function(epsilon,tau)/(1-htt)^(1/sqrt(tau*(1-tau))))
    }
  } else {

    ## Categorical predictors - this is the workhorse
    z <- as.matrix(z)
    num.z <- NCOL(z)
    epsilon <- numeric(length=n)
    htt <- numeric(length=n)
    ## At least one predictor for which degree > 0
    if(any(K[,1] > 0)) {
      P <- prod.spline(x=x,K=K,knots=knots,basis=basis)
      ## Check for 1 or fewer degrees of freedom
      if(NCOL(P) >= (n-1)) return(sqrt(.Machine$double.xmax))
      for(i in 1:nrow.z.unique) {
        zz <- ind == ind.vals[i]
        L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
        if(basis=="additive" || basis=="glp") {
          ## Additive spline regression models have an intercept in the lm()
          ## model (though not in the gsl.bs function)
          if(is.null(tau)) {
            model <- lm.wfit(cbind(1,P),y,L)
            ## Check for rank-deficient fit
            if(model$rank != (NCOL(P)+1)) return(sqrt(.Machine$double.xmax))
          } else {
            model <- rq.wfit(cbind(1,P),y,weights=L,tau=tau)
            ## Check for rank-deficient fit
            if(NCOL(model$x) != (NCOL(P)+1)) return(sqrt(.Machine$double.xmax))
          }
        } else {
          if(is.null(tau)) {
            model <- lm.wfit(P,y,L)
            ## Check for rank-deficient fit          
            if(model$rank != NCOL(P)) return(sqrt(.Machine$double.xmax))
          } else {
            model <- rq.wfit(P,y,weights=L,tau=tau)
            ## Check for rank-deficient fit          
            if(NCOL(model$x) != NCOL(P)) return(sqrt(.Machine$double.xmax))
          }
        }
        epsilon[zz] <- residuals(model)[zz]
        if(is.null(tau))
          htt[zz] <- hat(model$qr)[zz]
        else
          htt[zz] <- hat(model$x)[zz]  ## XXXXX this is absolutely wrong...
      }
    } else {
      ## No predictors for which degree > 0      
      z.factor <- data.frame(factor(z[,1]))
      if(num.z > 1) for(i in 2:num.z) z.factor <- data.frame(z.factor,factor(z[,i]))
      for(i in 1:nrow.z.unique) {
        zz <- ind == ind.vals[i]
        L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda,is.ordered.z=is.ordered.z)
        ## Whether we use additive, glp, or tensor products, this
        ## model has no continuous predictors hence the intercept is
        ## the parameter that may shift with the categorical
        ## predictors
        if(is.null(tau)) {
          model <- lm.wfit(matrix(1,n,1),y,L)
          htt[zz] <- hat(model$qr)[zz]
        } else {
          model <- rq.wfit(matrix(1,n,1),y,weights=L,tau=tau)
          htt[zz] <- hat(model$x)[zz] ## XXXXX This most certainly is incorrect! But it is a corner case so ignore for now 6/11/11 (all x irrelevant)
        }          
        epsilon[zz] <- residuals(model)[zz]
      }
    }
    
    htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
    
    if(cv.func == "cv.ls") {
      cv <- mean((epsilon/(1-htt))^2)
    } else if(cv.func == "cv.gcv"){
      cv <- mean((epsilon/(1-mean(htt)))^2)
    } else if(cv.func == "cv.aic"){
      sigmasq <- mean(epsilon^2)
      traceH <- sum(htt)
      penalty <- (1+traceH/n)/(1-(traceH+2)/n)
      cv <- ifelse(penalty < 0, .Machine$double.xmax, log(sigmasq)+penalty);
    } else if(cv.func == "cv.rq") {
      ## Note - this is defined in util.R so if you modify there you must modify here also
      cv <- mean(check.function(epsilon,tau)/(1-htt)^(1/sqrt(tau*(1-tau))))
    }
    
  }
  
  return(cv)
  
}

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
																		 cv.func=c("cv.ls","cv.gcv","cv.aic","cv.rq"),
                                     tau=NULL) {

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
                           tau=tau)

    cv.uniform <- cv.factor.spline(x=x, 
                                   y=y, 
                                   z=z, 
                                   K=K, 
                                   I=I, 
                                   knots="uniform", 
                                   basis=basis, 
                                   cv.func=cv.func,
                                   tau=tau)
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
                           tau=tau)
    
  }

  attr(cv, "knots.opt") <- knots.opt
  
  return(cv)
}

## Feb 1 2011 - replaced lm() and hatvalues() with lsfit and hat. On
## large datasets makes a noticeable difference (1 predictor, n=10,000
## drops from 17 secs to 13, n=100,000 drops 92 to to 71 seconds, so
## cuts runtime to ~ 0.77 of original time)

cv.factor.spline <- function(x,
														 y,
														 z=NULL,
														 K,
														 I=NULL,
														 knots=c("quantiles","uniform"),
														 basis=c("additive","tensor","glp"),
														 cv.func=c("cv.ls","cv.gcv","cv.aic","cv.rq"),
                             tau=NULL) {
  
  if(missing(x) || missing(y) || missing (K)) stop(" must provide x, y and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")  
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  
  cv.func <- match.arg(cv.func)
  
  n <- NROW(x)
  
  ## Otherwise, compute the cross-validation function
  
  if(any(K[,1] > 0)||any(I > 0)) {
    P <- prod.spline(x=x,z=z,K=K,I=I,knots=knots,basis=basis)
    ## Check for 1 or fewer degrees of freedom
    if(NCOL(P) >= (n-1)) return(sqrt(.Machine$double.xmax))
    if(basis=="additive" || basis=="glp") {
      ## Additive spline regression models have an intercept in the
      ## lm() model (though not in the gsl.bs function)
      if(is.null(tau)) {
        model <- lm.fit(cbind(1,P),y)
        ## Check for rank-deficient fit
        if(model$rank != (NCOL(P)+1)) return(sqrt(.Machine$double.xmax))
      } else {
        model <- rq.fit(cbind(1,P),y,tau=tau)
        ## Check for rank-deficient fit
        if(NCOL(model$x) != (NCOL(P)+1)) return(sqrt(.Machine$double.xmax))
      }
    } else {
      if(is.null(tau)) {
        model <- lm.fit(P,y)
        ## Check for rank-deficient fit      
        if(model$rank != NCOL(P)) return(sqrt(.Machine$double.xmax))
      } else {
        model <- rq.fit(P,y,tau=tau)
        ## Check for rank-deficient fit      
        if(NCOL(model$rx) != NCOL(P)) return(sqrt(.Machine$double.xmax))
      }
    }
    if(is.null(tau))
      htt <- hat(model$qr)
    else
      htt <- hat(model$x) ### XXX This most certainly is wrong! Quantiles may need further refinement...
    htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
    epsilon <- residuals(model)
  } else {
    htt <- rep(1/n,n)
    epsilon <- y-mean(y)
  }
  
  if(cv.func == "cv.ls") {
    cv <- mean((epsilon/(1-htt))^2)
  } else if(cv.func == "cv.gcv"){
    cv <- mean((epsilon/(1-mean(htt)))^2)
  } else if(cv.func == "cv.aic"){
    sigmasq <- mean(epsilon^2)
    traceH <- sum(htt)
    penalty <- (1+traceH/n)/(1-(traceH+2)/n)
    cv <- ifelse(penalty < 0, .Machine$double.xmax, log(sigmasq)+penalty);
  } else if(cv.func == "cv.rq") {
    ## Note - this is defined in util.R so if you modify there you must modify here also
    cv <- mean(check.function(epsilon,tau)/(1-htt)^(1/sqrt(tau*(1-tau))))
  }
  
  return(cv)
  
}
