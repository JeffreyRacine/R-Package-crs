## This file contains code for multivariate generalized local
## polynomial kernel regression with mixed datatypes. It relies on
## functions in the np package and on snomadr which currently resides
## in the crs package (June 29 2011).

## Note that the approach taken here is computationally efficient and
## relies on expressing the local polynomial method in slightly
## different form purely for computational simplicity. Both approaches
## are identical though.

NZD <- function(a) {
  sapply(1:NROW(a), function(i) {if(a[i] < 0) min(-.Machine$double.xmin,a[i]) else max(.Machine$double.xmin,a[i])})
}

mypoly <- function(x,degree,gradient.compute=FALSE,r=0) {

  if(missing(x)) stop(" Error: x required")
  if(missing(degree)) stop(" Error: degree required")
  if(degree < 1) stop(" Error: degree must be a positive integer")
  n <- degree + 1

  if(gradient.compute) {
    Z <- NULL
    for(i in 1:degree) {
      if((i-r) >= 0) {
        tmp <- (factorial(i)/factorial(i-r))*x^max(0,i-r)
      } else {
        tmp <- rep(0,length(x))
      }
      Z <- cbind(Z,tmp)
    }
    ## If we want to send back matrix of zeros (i.e will not count in
    ## computation of derivative) return this baby.
    if(r == -1) Z <- matrix(0,NROW(Z),NCOL(Z))
  } else {
      Z <- outer(x,1L:degree,"^")
  }
  
  return(as.matrix(Z))

}

## W.glp is a modified version of the polym() function (stats). The
## function accepts a vector of degrees and provides a generalized
## polynomial with varying polynomial order.

W.glp <- function(xdat = NULL,
                  degree = NULL,
                  gradient.vec = NULL) {

  if(is.null(xdat)) stop(" Error: You must provide data")
  if(is.null(degree) | any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  xdat <- as.data.frame(xdat)

  xdat.col.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})
  k <- ncol(as.data.frame(xdat[,xdat.col.numeric]))

  if(k > 0) {
    xdat.numeric <- as.data.frame(xdat[,xdat.col.numeric])
  }

  num.numeric <- ncol(xdat.numeric)

  if(!is.null(gradient.vec) && (length(gradient.vec) != num.numeric)) stop(paste(" Error: gradient vector and number of numeric predictors must be conformable\n",sep=""))
  if(!is.null(gradient.vec) && any(gradient.vec < 0)) stop(paste(" Error: gradient vector must contain non-negative integers\n",sep=""))
  if(!is.null(gradient.vec)) gradient.vec[gradient.vec==0] <- -1
  if(!is.null(gradient.vec)) {
    gradient.compute <- TRUE
  } else {
    gradient.compute <- FALSE
    gradient.vec <- rep(NA,NCOL(xdat))
  }

  if(length(degree) != num.numeric) stop(" Error: degree vector and number of numeric predictors incompatible")

  if(all(degree == 0) | k == 0) {

    ## Local constant OR no continuous variables

    return(matrix(1,nrow=nrow(xdat.numeric),ncol=1))

  } else {

    degree.list <- list()
    for(i in 1:k) degree.list[[i]] <- 0:degree[i]
    z <- do.call("expand.grid", degree.list, k)
    s <- rowSums(z)
    ind <- (s > 0) & (s <= max(degree))
    z <- z[ind, ,drop=FALSE]
    if(!all(degree==max(degree))) {
      for(j in 1:length(degree)) {
        d <- degree[j]
        if((d < max(degree)) & (d > 0)) {
          s <- rowSums(z)
          d <- (s > d) & (z[,j,drop=FALSE]==matrix(d,nrow(z),1,byrow=TRUE))
          z <- z[!d, ]
        }
      }
    }
    res <- rep.int(1,nrow(xdat.numeric))
    if(degree[1] > 0) res <- cbind(1, mypoly(xdat.numeric[,1],degree[1],gradient.compute=gradient.compute,r=gradient.vec[1]))[, 1 + z[, 1]]
    if(k > 1) for (i in 2:k) if(degree[i] > 0) res <- res * cbind(1, mypoly(xdat.numeric[,i],degree[i],gradient.compute=gradient.compute,r=gradient.vec[i]))[, 1 + z[, i]]
    res <- matrix(res,nrow=NROW(xdat))
    colnames(res) <- apply(z, 1L, function(x) paste(x, collapse = "."))
    return(as.matrix(cbind(1,res)))

  }

}

npglpreg <- function(...) UseMethod("npglpreg")

npglpreg.default <- function(tydat=NULL,
                             txdat=NULL,
                             eydat=NULL,
                             exdat=NULL,
                             bws=NULL,
                             degree=NULL,
                             leave.one.out=FALSE,
                             ukertype=c("liracine","aitchisonaitken"),
                             okertype=c("liracine","wangvanryzin"),
                             bwtype = c("fixed","generalized_nn","adaptive_nn"),
                             gradient.vec=NULL,
                             gradient.categorical=FALSE,
                             ...) {

  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)

  est <- glpregEst(tydat=tydat,
                   txdat=txdat,
                   eydat=eydat,
                   exdat=exdat,
                   bws=bws,
                   degree=degree,
                   leave.one.out=leave.one.out,
                   ukertype=ukertype,
                   okertype=okertype,
                   bwtype=bwtype,
                   gradient.vec=gradient.vec,
                   ...)

  ## Gradients for categorical predictors (here we just compute them
  ## all though gradient for continuous predictors given above
  ## requires specific variables to be selected)

  est$gradient.categorical.mat <- NULL

  if(gradient.categorical) {

    num.numeric <- est$num.numeric
    num.categorical <- est$num.categorical

    if(num.categorical > 0) {

      num.eval <- ifelse(is.null(exdat),nrow(txdat),nrow(exdat))
      gradient.categorical.mat <- matrix(NA,nrow=num.eval,ncol=num.categorical)

      for(i in 1:num.categorical) {

        if(is.null(exdat)) {
          exdat.base <- txdat
        } else {
          exdat.base <- exdat
        }

        categorical.index <- est$categorical.index
        eval.base <- levels(txdat[,categorical.index[i]])[1]
        eval.levels <- levels(txdat[,categorical.index[i]])

        if(is.ordered(txdat[,categorical.index[i]])) {
          exdat.base[,categorical.index[i]] <- ordered(rep(eval.base,num.eval),levels=eval.levels)
        } else {
          exdat.base[,categorical.index[i]] <- factor(rep(eval.base,num.eval),levels=eval.levels)
        }

        est.base <- glpregEst(tydat=tydat,
                              txdat=txdat,
                              eydat=eydat,
                              exdat=exdat.base,
                              bws=bws,
                              degree=degree,
                              leave.one.out=leave.one.out,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,
                              gradient.vec=NULL,
                              ...)

        gradient.categorical.mat[,i] <- est$fitted.values - est.base$fitted.values

      }

      est$gradient.categorical.mat <- gradient.categorical.mat

    }

  }

  ## Add results to estimated object.

  if(!is.null(eydat)) {
    est$r.squared <- RSQfunc(eydat,est$fitted.values)
    est$residuals <- eydat - est$fitted.values
  } else if(is.null(eydat)&&is.null(exdat)) {
    est$r.squared <- RSQfunc(tydat,est$fitted.values)
    est$residuals <- tydat - est$fitted.values
  } else {
    est$r.squared <- NULL
    est$residuals <- NULL
  }
    
  est$call <- match.call()

  ## Return object of type npglpreg

  class(est) <- "npglpreg"

  return(est)

}

## Basic print method.

print.npglpreg <- function(x,
                           ...) {
  cat("Call:\n")
  print(x$call)

}

summary.npglpreg <- function(object,
                             ...) {

  cat("Call:\n")
  print(object$call)
  cat("\nGeneralized Local Polynomial Kernel Regression\n",sep="")

  ## Summarize continuous predictors

  if(object$num.numeric == 1){
    cat(paste("\nThere is ",format(object$num.numeric), " continuous predictor",sep=""),sep="")
  } else if(object$num.numeric > 1) {
    cat(paste("\nThere are ",format(object$num.numeric), " continuous predictors",sep=""),sep="")
  }

  if(object$num.numeric > 1) for(j in 1:object$num.numeric) 
      cat(paste("\nBandwidth for ",format(object$xnames[object$numeric.index][j]),": ",format(object$bws[object$numeric.index][j]),sep=""),sep="")

  for(j in 1:object$num.numeric)
    cat(paste("\nDegree for ",format(object$xnames[object$numeric.index][j]),": ",format(object$degree[j]),sep=""),sep="")

  ## Summarize categorical predictors  
    
  if(object$num.categorical==1) {
    cat(paste("\nThere is ",format(object$num.categorical), " categorical predictor",sep=""),sep="")
  } else if(object$num.categorical > 1) {
    cat(paste("\nThere are ",format(object$num.categorical), " categorical predictors",sep=""),sep="")
  }

  if(object$num.categorical > 1) for(j in 1:(object$num.categorical)) 
    cat(paste("\nBandwidth for ",format(object$xnames[object$categorical.index][j]),": ",format(object$bws[object$categorical.index][j]),sep=""),sep="")


  ## Summary statistics

  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4), sep=""))
  if(!is.null(object$fv)) {
    cat(paste("\nCross-validation score: ", format(object$fv,digits=8), sep=""))
    cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
  }
  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))
  cat("\n\n")

}

## Method for predicting given a new data frame.

predict.npglpreg <- function(object,
                             newdata=NULL,
                             gradient.vec=NULL,
                             ...) {

  if(!is.null(newdata)&&nrow(newdata)==1) stop(" Error: newdata must have more than one row")
  
  if(is.null(newdata)) {
    
    ## If no new data provided, return sample fit.
    fitted.values <- fitted(object)
    gradient <- object$gradient
    gradient.categorical.mat <- object$gradient.categorical.mat
    
  } else{

    ## Get training data from object (xz and y) and parse into factors
    ## and numeric.
    
    degree <- object$degree
    bws <- object$bws
    bwtype <- object$bwtype
    ukertype <- object$ukertype
    okertype <- object$okertype
    if(is.null(gradient.vec)) {
      gradient.vec <- object$gradient.vec
    } 
    
    txdat <- object$x
    tydat <- object$y

    tt <- terms(object)
    has.ey <- succeedWithResponse(tt, newdata)
    if (has.ey) {
      eydat <- model.response(model.frame(tt,newdata))
    } else {
      eydat <- NULL
    }      
    exdat <- model.frame(delete.response(tt),newdata,xlev=object$xlevels)

    ## Return the predicted values.

    est <- npglpreg.default(tydat=tydat,
                            txdat=txdat,
                            exdat=exdat,
                            eydat=eydat,
                            bws=bws,
                            degree=degree,
                            ukertype=ukertype,
                            okertype=okertype,
                            bwtype=bwtype,
                            gradient.vec=gradient.vec,
                            ...)
    
    fitted.values <- est$fitted.values
    gradient <- est$gradient
    gradient.categorical.mat <- est$gradient.categorical.mat
    
  }

  attr(fitted.values, "gradient") <- gradient
  attr(fitted.values, "gradient.categorical.mat") <- gradient.categorical.mat

  return(fitted.values)

}

npglpreg.formula <- function(formula,
                             data=list(),
                             tydat=NULL,
                             txdat=NULL,
                             eydat=NULL,
                             exdat=NULL,
                             bws=NULL,
                             degree=NULL,
                             leave.one.out=FALSE,
                             ukertype=c("liracine","aitchisonaitken"),
                             okertype=c("liracine","wangvanryzin"),
                             bwtype = c("fixed","generalized_nn","adaptive_nn"),
                             cv=c("degree-bandwidth","bandwidth","none"),
                             cv.func=c("cv.ls","cv.gcv","cv.aic"),
                             opts=list("MAX_BB_EVAL"=10000,
                               "EPSILON"=.Machine$double.eps,
                               "INITIAL_MESH_SIZE"="1.0e-01",
                               "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                               "MIN_POLL_SIZE"=paste("r",sqrt(.Machine$double.eps),sep="")),
                             nmulti=5,
                             random.seed=42,
                             degree.max=5,
                             degree.min=0,
                             bandwidth.max=1.0e+05,
                             bandwidth.min=1.0e-03,
                             gradient.vec=NULL,
                             gradient.categorical=FALSE,
                             ridge.warning=FALSE,
                             ...) {

  if(!require(np)) stop(" Error: you must install the np package to use this function")

  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)
  cv <- match.arg(cv)
  cv.func <- match.arg(cv.func)

  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  tydat <- model.response(mf)
  txdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

  fv <- NULL
  ptm <- system.time("")

  if(cv!="none") {
    ptm <- ptm + system.time(model.cv <-glpcvNOMAD(ydat=tydat,
                                                   xdat=txdat,
                                                   opts=opts,
                                                   cv=cv,
                                                   degree=degree,
                                                   bandwidth=bws,
                                                   bwmethod=cv.func,
                                                   bwtype=bwtype,
                                                   nmulti=nmulti,
                                                   random.seed=random.seed,
                                                   degree.max=degree.max,
                                                   degree.min=degree.min,
                                                   bandwidth.max=bandwidth.max,
                                                   bandwidth.min=bandwidth.min,
                                                   ...))
    degree <- model.cv$degree
    bws <- model.cv$bws
    fv <- model.cv$fv
  }
  
  ptm <- ptm + system.time(est <- npglpreg.default(tydat=tydat,
                                                   txdat=txdat,
                                                   eydat=eydat,
                                                   exdat=exdat,
                                                   bws=bws,
                                                   degree=degree,
                                                   leave.one.out=leave.one.out,
                                                   ukertype=ukertype,
                                                   okertype=okertype,
                                                   bwtype=bwtype,
                                                   gradient.vec=gradient.vec,
                                                   gradient.categorical=gradient.categorical,
                                                   ...))
  
  est$call <- match.call()
  est$formula <- formula
  est$terms <- mt
  est$xlevels <- .getXlevels(mt, mf)
  est$x <- txdat
  est$y <- tydat
  est$fv <- fv
  est$nmulti <- nmulti
  est$ptm <- ptm

  return(est)

}

glpregEst <- function(tydat=NULL,
                      txdat=NULL,
                      eydat=NULL,
                      exdat=NULL,
                      bws=NULL,
                      degree=NULL,
                      leave.one.out=FALSE,
                      ukertype=c("liracine","aitchisonaitken"),
                      okertype=c("liracine","wangvanryzin"),
                      bwtype=c("fixed","generalized_nn","adaptive_nn"),
                      gradient.vec=NULL,
                      ridge.warning=FALSE,
                      ...) {

  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)

  if(is.null(tydat)) stop(" Error: You must provide y data")
  if(is.null(txdat)) stop(" Error: You must provide X data")
  if(is.null(bws)) stop(" Error: You must provide a bandwidth object")
  if(is.null(degree) | any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  miss.ex = is.null(exdat)
  miss.ey = is.null(eydat)

  if (miss.ex){
    exdat <- txdat
  }

  txdat <- as.data.frame(txdat)
  exdat <- as.data.frame(exdat)

  maxPenalty <- sqrt(.Machine$double.xmax)

  n.train <- nrow(txdat)
  n.eval <- nrow(exdat)

  xdat.numeric <- sapply(1:ncol(txdat),function(i){is.numeric(txdat[,i])})
  categorical.index <- which(xdat.numeric==FALSE)
  numeric.index <- which(xdat.numeric==TRUE)  
  num.numeric <- sum(sapply(1:NCOL(txdat),function(i){is.numeric(txdat[,i])})==TRUE)
  num.categorical <- NCOL(txdat)-num.numeric

  ## Check whether it appears that training and evaluation data are
  ## conformable

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  if(ncol(txdat)!=ncol(exdat))
    stop(" Error: training and evaluation data have unequal number of columns\n")

  if(all(degree == 0)) {

    ## Local constant using only one call to npksum

    if(leave.one.out == TRUE) {

      ## exdat not supported with leave.one.out, but this is only used
      ## for cross-validation hence no exdat

      tww <- npksum(txdat = txdat,
                    weights = as.matrix(data.frame(1,tydat)),
                    tydat = rep(1,length(tydat)),
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

    } else {

      tww <- npksum(txdat = txdat,
                    exdat = exdat,
                    weights = as.matrix(data.frame(1,tydat)),
                    tydat = rep(1,length(tydat)),
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

    }


    ## Note that as bandwidth approaches zero the local constant
    ## estimator undersmooths and approaches each sample realization,
    ## so use the convention that when the sum of the kernel weights
    ## equals 0, return y. This is unique to this code.

    mhat <- tww[2,]/NZD(tww[1,])

    return(list(fitted.values = mhat,
                gradient = NULL,
                coef.mat = NULL,
                bwtype = bwtype,
                ukertype = ukertype,
                okertype = okertype,
                degree = degree,
                bws = bws,
                nobs = n.train,
                num.numeric = num.numeric,
                num.categorical = num.categorical,
                xnames = names(txdat),
                categorical.index = categorical.index,
                numeric.index = numeric.index,
                gradient.vec = gradient.vec))
    
  } else {

    W <- W.glp(txdat,degree)

    if(miss.ex) {
      W.eval <- W
    } else {
      W.eval <- W.glp(exdat,degree)
    }
    
    if(!is.null(gradient.vec)) W.eval.deriv <- W.glp(exdat,degree,gradient.vec=gradient.vec)

    ## Local polynomial via smooth coefficient formulation and one
    ## call to npksum

    if(leave.one.out == TRUE) {

      ## exdat not supported with leave.one.out, but this is only used
      ## for cross-validation hence no exdat

      tww <- npksum(txdat = txdat,
                    tydat = as.matrix(cbind(tydat,W)),
                    weights = W,
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

    } else {

      tww <- npksum(txdat = txdat,
                    exdat = exdat,
                    tydat = as.matrix(cbind(tydat,W)),
                    weights = W,
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

    }

    tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n.eval))[1,,]

    ## June 30 2011... spent hours trying to track down an issue with
    ## this code. When there is only one row in the evaluation data
    ## this array becomes a matrix and things break. Should be simple
    ## to fix but not so simple it appears... this needs attention for
    ## a robust package, but for now I trap this in predict.glpregEst
    ## and throw an error to avoid issues in the future...
    tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n.eval))[-1,,]

    coef.mat <- matrix(maxPenalty,ncol(W),n.eval)
    epsilon <- 1.0/n.eval
    ridge <- double(n.eval)
    doridge <- !logical(n.eval)

    nc <- ncol(tww[,,1])

    ridger <- function(i) {
      doridge[i] <<- FALSE
      ridge.val <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
      tryCatch(solve(tww[,,i]+diag(rep(ridge[i],nc)),
                     tyw[,i]+c(ridge.val,rep(0,nc-1))),
               error = function(e){
                 ridge[i] <<- ridge[i]+epsilon
                 doridge[i] <<- TRUE
                 if(ridge.warning)
                   console <- printPush(paste("\rWarning: ridging required for inversion at obs. = ", i, ", ridge = ",formatC(ridge[i],digits=4,format="f"),"        ",sep=""),console = console)
                 return(rep(maxPenalty,nc))
               })
    }

    ## Test for singularity of the generalized local polynomial
    ## estimator, shrink the mean towards the local constant mean.

    while(any(doridge)){
      iloo <- (1:n.eval)[doridge]
      coef.mat[,iloo] <- sapply(iloo, ridger)
    }

    mhat <- sapply(1:n.eval, function(i) {
      W.eval[i,, drop = FALSE] %*% coef.mat[,i]
    })

    if(!is.null(gradient.vec)) {
      gradient <- sapply(1:n.eval, function(i) {
        W.eval.deriv[i,-1, drop = FALSE] %*% coef.mat[-1,i]
      })
    } else {
      gradient <- NULL
    }

    return(list(fitted.values = mhat,
                gradient = gradient,
                coef.mat = t(coef.mat[-1,]),
                bwtype = bwtype,
                ukertype = ukertype,
                okertype = okertype,
                degree = degree,
                bws = bws,
                nobs = n.train,
                num.numeric = num.numeric,
                num.categorical = num.categorical,
                xnames = names(txdat),
                categorical.index = categorical.index,
                numeric.index = numeric.index,
                gradient.vec = gradient.vec))
    
  }

}

minimand.cv.ls <- function(bws=NULL,
                           ydat=NULL,
                           xdat=NULL,
                           degree=NULL,
                           W=NULL,
                           ukertype=c("liracine","aitchisonaitken"),
                           okertype=c("liracine","wangvanryzin"),
                           bwtype = c("fixed","generalized_nn","adaptive_nn"),
                           ridge.warning=FALSE,
                           ...) {

  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)  

  if(is.null(ydat)) stop(" Error: You must provide y data")
  if(is.null(xdat)) stop(" Error: You must provide X data")
  if(is.null(W)) stop(" Error: You must provide a weighting matrix W")
  if(is.null(bws)) stop(" Error: You must provide a bandwidth object")
  if(is.null(degree) | any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  xdat <- as.data.frame(xdat)

  n <- length(ydat)

  maxPenalty <- sqrt(.Machine$double.xmax)

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  if(any(bws<=0)) {

    return(maxPenalty)

  } else {

    if(all(degree == 0)) {

      ## Local constant via one call to npksum

      tww <- npksum(txdat = xdat,
                    weights = as.matrix(data.frame(1,ydat)),
                    tydat = rep(1,n),
                    bws = bws,
                    leave.one.out = TRUE,
                    bandwidth.divide = TRUE,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

      mean.loo <- tww[2,]/NZD(tww[1,])
      
      if (!any(mean.loo == maxPenalty)){
        fv <- mean((ydat-mean.loo)^2)
      } else {
        fv <- maxPenalty
      }

      fv <- ifelse(is.finite(fv),fv,maxPenalty)

      console <- printPush("\r                                                ",console = console)
      console <- printPush(paste("\rfv = ",format(fv)," ",sep=""),console = console)

      return(fv)

    } else {


      ## Generalized local polynomial via smooth coefficient
      ## formulation and one call to npksum

      tww <- npksum(txdat = xdat,
                    tydat = as.matrix(cbind(ydat,W)),
                    weights = W,
                    bws = bws,
                    leave.one.out = TRUE,
                    bandwidth.divide = TRUE,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

      tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[1,,]
      tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[-1,,]

      mean.loo <- rep(maxPenalty,n)
      epsilon <- 1.0/n
      ridge <- double(n)
      doridge <- !logical(n)

      nc <- ncol(tww[,,1])

      ## Test for singularity of the generalized local polynomial
      ## estimator, shrink the mean towards the local constant mean.

      ridger <- function(i) {
        doridge[i] <<- FALSE
        ridge.val <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
        W[i,, drop = FALSE] %*% tryCatch(solve(tww[,,i]+diag(rep(ridge[i],nc)),
                tyw[,i]+c(ridge.val,rep(0,nc-1))),
                error = function(e){
                  ridge[i] <<- ridge[i]+epsilon
                  doridge[i] <<- TRUE
                  if(ridge.warning)
                    console <- printPush(paste("\rWarning: ridging required for inversion at obs. = ", i, ", ridge = ",formatC(ridge[i],digits=4,format="f"),"        ",sep=""),console = console)
                  return(rep(maxPenalty,nc))
                })
      }

      while(any(doridge)){
        iloo <- (1:n)[doridge]
        mean.loo[iloo] <- sapply(iloo, ridger)
      }

      if (!any(mean.loo == maxPenalty)){
        fv <- mean((ydat-mean.loo)^2)
      } else {
        fv <- maxPenalty
      }

      fv <- ifelse(is.finite(fv),fv,maxPenalty)

      console <- printPush("\r                                                ",console = console)
      console <- printPush(paste("\rfv = ",format(fv)," ",sep=""),console = console)

      return(fv)

    }

  }

}

minimand.cv.aic <- function(bws=NULL,
                            ydat=NULL,
                            xdat=NULL,
                            degree=NULL,
                            W=NULL,
                            ukertype=c("liracine","aitchisonaitken"),
                            okertype=c("liracine","wangvanryzin"),
                            bwtype = c("fixed","generalized_nn","adaptive_nn"),
                            ridge.warning=FALSE,
                            ...) {

  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)  

  if(is.null(ydat)) stop(" Error: You must provide y data")
  if(is.null(xdat)) stop(" Error: You must provide X data")
  if(!all(degree==0)) if(is.null(W)) stop(" Error: You must provide a weighting matrix W")
  if(is.null(bws)) stop(" Error: You must provide a bandwidth object")
  if(is.null(degree) | any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  xdat <- as.data.frame(xdat)

  n <- length(ydat)

  maxPenalty <- sqrt(.Machine$double.xmax)

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  if(any(bws<=0)) {

    return(maxPenalty)

  } else {

    ## This computes the kernel function when i=j (i.e., K(0))

    kernel.i.eq.j <- npksum(txdat = xdat[1,],
                            weights = as.matrix(data.frame(1,ydat)[1,]),
                            tydat = 1,
                            bws = bws,
                            bandwidth.divide = TRUE,
                            ukertype = ukertype,
                            okertype = okertype,
                            bwtype = bwtype,                    
                            ...)$ksum[1,1]

    if(all(degree == 0)) {

      ## Local constant via one call to npksum

      tww <- npksum(txdat = xdat,
                    weights = as.matrix(data.frame(1,ydat)),
                    tydat = rep(1,n),
                    bws = bws,
                    bandwidth.divide = TRUE,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

      ghat <- tww[2,]/NZD(tww[1,])

      trH <- kernel.i.eq.j*sum(1/NZD(tww[1,]))

      aic.penalty <- (1+trH/n)/(1-(trH+2)/n)

      if (!any(ghat == maxPenalty) & (aic.penalty > 0)){
        fv <- log(mean((ydat-ghat)^2)) + aic.penalty
      } else {
        fv <- maxPenalty
      }

      return(ifelse(is.finite(fv),fv,maxPenalty))

    } else {

      ## Generalized local polynomial via smooth coefficient
      ## formulation and one call to npksum

      tww <- npksum(txdat = xdat,
                    tydat = as.matrix(cbind(ydat,W)),
                    weights = W,
                    bws = bws,
                    bandwidth.divide = TRUE,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,                    
                    ...)$ksum

      tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[1,,]
      tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[-1,,]

      ghat <- rep(maxPenalty,n)
      epsilon <- 1.0/n
      ridge <- double(n)
      doridge <- !logical(n)

      nc <- ncol(tww[,,1])

      ## Test for singularity of the generalized local polynomial
      ## estimator, shrink the mean towards the local constant mean.

      ridger <- function(i) {
        doridge[i] <<- FALSE
        ridge.val <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
        W[i,, drop = FALSE] %*% tryCatch(solve(tww[,,i]+diag(rep(ridge[i],nc)),
                tyw[,i]+c(ridge.val,rep(0,nc-1))),
                error = function(e){
                  ridge[i] <<- ridge[i]+epsilon
                  doridge[i] <<- TRUE
                  if(ridge.warning)
                    console <- printPush(paste("\rWarning: ridging required for inversion at obs. = ", i, ", ridge = ",formatC(ridge[i],digits=4,format="f"),"        ",sep=""),console = console)
                  return(rep(maxPenalty,nc))
                })
      }

      while(any(doridge)){
        ii <- (1:n)[doridge]
        ghat[ii] <- sapply(ii, ridger)
      }

      trH <- kernel.i.eq.j*sum(sapply(1:n,function(i){
        W[i,, drop = FALSE] %*% solve(tww[,,i]+diag(rep(ridge[i],nc))) %*% t(W[i,, drop = FALSE])
      }))

      if (!any(ghat == maxPenalty)){
        fv <- log(mean((ydat-ghat)^2)) + (1+trH/n)/(1-(trH+2)/n)
      } else {
        fv <- maxPenalty
      }

      fv <- ifelse(is.finite(fv),fv,maxPenalty)

      console <- printPush("\r                                                ",console = console)
      console <- printPush(paste("\rfv = ",format(fv)," ",sep=""),console = console)

      return(fv)

    }

  }

}

glpcv <- function(ydat=NULL,
                  xdat=NULL,
                  degree=NULL,
                  bwmethod=c("cv.ls","cv.aic"),
                  ukertype=c("liracine","aitchisonaitken"),
                  okertype=c("liracine","wangvanryzin"),
                  bwtype = c("fixed","generalized_nn","adaptive_nn"),
                  nmulti=NULL,
                  random.seed=42,
                  optim.maxattempts = 10,
                  optim.method=c("Nelder-Mead", "BFGS", "CG"),
                  optim.reltol=sqrt(.Machine$double.eps),
                  optim.abstol=.Machine$double.eps,
                  optim.maxit=500,
                  debug=FALSE,
                  ...) {

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  if(debug) system("rm optim.debug bandwidth.out optim.out")

  if(is.null(ydat)) stop(" Error: You must provide y data")
  if(is.null(xdat)) stop(" Error: You must provide X data")
  if(is.null(degree) | any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))
  if(!is.null(nmulti) && nmulti < 1) stop(paste(" Error: nmulti must be a positive integer (minimum 1)\nnmulti is (", nmulti, ")\n",sep=""))

  bwmethod <- match.arg(bwmethod)

  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)  

  optim.method <- match.arg(optim.method)
  optim.control <- list(abstol = optim.abstol,
                        reltol = optim.reltol,
                        maxit = optim.maxit)

  maxPenalty <- sqrt(.Machine$double.xmax)

  xdat <- as.data.frame(xdat)

  num.bw <- ncol(xdat)

  if(is.null(nmulti)) nmulti <- min(5,num.bw)

  ## Which variables are categorical, which are discrete...

  xdat.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})

  ## First initialize initial search values of the vector of
  ## bandwidths to lie in [0,1]

  if(debug) write(c("cv",paste(rep("x",num.bw),seq(1:num.bw),sep="")),file="optim.debug",ncol=(num.bw+1))

  ## Pass in the local polynomial weight matrix rather than
  ## recomputing with each iteration.

  W <- W.glp(xdat,degree)

  sum.lscv <- function(bw.gamma,...) {

    ## Note - we set the kernel for unordered and ordered regressors
    ## to the liracine kernel (0<=lambda<=1) and test for proper
    ## bounds in sum.lscv.

    if(all(bw.gamma>=0)&&all(bw.gamma[!xdat.numeric]<=1)) {
      lscv <- minimand.cv.ls(bws=bw.gamma,
                             ydat=ydat,
                             xdat=xdat,
                             ukertype=ukertype,
                             okertype=okertype,
                             bwtype=bwtype,                    
                             ...)
    } else {
      lscv <- maxPenalty
    }

    if(debug) write(c(lscv,bw.gamma),file="optim.debug",ncol=(num.bw+1),append=TRUE)
    return(lscv)
  }

  sum.aicc <- function(bw.gamma,...) {

    ## Note - we set the kernel for unordered and ordered regressors
    ## to the liracine kernel (0<=lambda<=1) and test for proper
    ## bounds in sum.lscv.

    if(all(bw.gamma>=0)&&all(bw.gamma[!xdat.numeric]<=1)) {
      aicc <- minimand.cv.aic(bws=bw.gamma,
                              ydat=ydat,
                              xdat=xdat,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,                    
                              ...)
    } else {
      aicc <- maxPenalty
    }

    if(debug) write(c(aicc,bw.gamma),file="optim.debug",ncol=(num.bw+1),append=TRUE)
    return(aicc)
  }

  ## Multistarting

  fv.vec <- numeric(nmulti)

  ## Pass in the W matrix rather than recomputing it each time

  for(iMulti in 1:nmulti) {

    num.numeric <- ncol(as.data.frame(xdat[,xdat.numeric]))

    ## First initialize to values for factors (`liracine' kernel)

    init.search.vals <- runif(ncol(xdat),0,1)

    for(i in 1:ncol(xdat)) {
      if(xdat.numeric[i]==TRUE) {
        init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
      }
    }

    ## Initialize `best' values prior to search

    if(iMulti == 1) {
      fv <- maxPenalty
      numimp <- 0
      bw.opt <- init.search.vals
      best <- 1
    }

    if(bwmethod == "cv.ls" ) {

      suppressWarnings(optim.return <- optim(init.search.vals,
                                             fn=sum.lscv,
                                             method=optim.method,
                                             control=optim.control,
                                             degree=degree,
                                             W=W,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,                    
                                             ...))

      attempts <- 0
      while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
        init.search.vals <- runif(ncol(xdat),0,1)
        if(xdat.numeric[i]==TRUE) {
          init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
        }
        attempts <- attempts + 1
        optim.control$abstol <- optim.control$abstol * 10.0
        optim.control$reltol <- optim.control$reltol * 10.0
        suppressWarnings(optim.return <- optim(init.search.vals,
                                               fn=sum.lscv,
                                               method=optim.method,
                                               control=optim.control,
                                               degree=degree,
                                               W=W,
                                               ukertype=ukertype,
                                               okertype=okertype,
                                               bwtype=bwtype,                    
                                               ...))
      }

    } else {

      suppressWarnings(optim.return <- optim(init.search.vals,
                                             fn=sum.aicc,
                                             method=optim.method,
                                             control=optim.control,
                                             degree=degree,
                                             W=W,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,                    
                                             ...))

      attempts <- 0
      while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
        init.search.vals <- runif(ncol(xdat),0,1)
        if(xdat.numeric[i]==TRUE) {
          init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
        }
        attempts <- attempts + 1
        optim.control$abstol <- optim.control$abstol * 10.0
        optim.control$reltol <- optim.control$reltol * 10.0
        suppressWarnings(optim.return <- optim(init.search.vals,
                                               fn = sum.aicc,
                                               method=optim.method,
                                               control = optim.control,
                                               W=W,
                                               ukertype=ukertype,
                                               okertype=okertype,
                                               bwtype=bwtype,                    
                                               ...))
      }
    }

    if(optim.return$convergence != 0) warning(" optim failed to converge")

    fv.vec[iMulti] <- optim.return$value

    if(optim.return$value < fv) {
      bw.opt <- optim.return$par
      fv <- optim.return$value
      numimp <- numimp + 1
      best <- iMulti
      if(debug) {
        if(iMulti==1) {
          write(cbind(iMulti,t(bw.opt)),"bandwidth.out",ncol=(1+length(bw.opt)))
          write(cbind(iMulti,fv),"optim.out",ncol=2)
        } else {
          write(cbind(iMulti,t(bw.opt)),"bandwidth.out",ncol=(1+length(bw.opt)),append=TRUE)
          write(cbind(iMulti,fv),"optim.out",ncol=2,append=TRUE)
        }
      }
    }

  }

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  return(list(bws=bw.opt,fv=fv,numimp=numimp,best=best,fv.vec=fv.vec))

}

## Note bws differ dramatically for numeric/categorical predictors and
## lower bounds need to also be different... can be 0 for categorical,
## but not for numeric...

## Note for scaling used for initial bandwidths may need to adjust
## bandwidth.min (now fraction of sd) by n^{-1/(2p+1)} etc.

## Note - numeric.scale=TRUE rescales numeric predictors (not a robust
## rescaling which we could achieve via sweep

## If the degree or bandwidth are fed in they are used as initial values
## for search where appropriate.

## set bandwidth.max to be ten thousand times the standard deviation
## set bandwidth.min minimum bw is 1/10 standard deviation

glpcvNOMAD <- function(ydat=NULL,
                       xdat=NULL,
                       degree=NULL,
                       bandwidth=NULL,
                       bwmethod=c("cv.ls","cv.aic"),
                       ukertype=c("liracine","aitchisonaitken"),
                       okertype=c("liracine","wangvanryzin"),
                       bwtype = c("fixed","generalized_nn","adaptive_nn"),
                       cv=c("degree-bandwidth", "bandwidth"),
                       nmulti=NULL,
                       numeric.scale=TRUE,
                       random.seed=42,
                       degree.max=5,
                       degree.min=0,
                       bandwidth.max=1.0e+05,
                       bandwidth.min=1.0e-03,
                       opts=list(),
                       ...) {

  ## Save the seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)  
  bwmethod <- match.arg(bwmethod)
  cv <- match.arg(cv)

  if(is.null(ydat)) stop(" Error: You must provide y data")
  if(is.null(xdat)) stop(" Error: You must provide X data")
  if(!is.null(nmulti) && nmulti < 1) stop(paste(" Error: nmulti must be a positive integer (minimum 1)\nnmulti is (", nmulti, ")\n",sep=""))

  if(!is.null(degree) && any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  if(is.null(bandwidth.max)) bandwidth.max <- .Machine$double.xmax

  if(degree.min < 0 ) stop(" Error: degree.min must be a non-negative integer")
  if(degree.max < degree.min) stop(" Error: degree.max must exceed degree.min")

  if(bandwidth.min < 0) stop(" Error: bandwidth.min must be non-negative")
  if(bandwidth.max < bandwidth.min) stop(" Error: bandwidth.max must exceed bandwidth.min")

  if(cv=="degree-bandwidth") {
    if(!is.null(degree) && any(degree>degree.max)) stop(" Error: degree supplied but exceeds degree.max")
    if(!is.null(degree) && any(degree<degree.min)) stop(" Error: degree supplied but less than degree.min")
  }
  
  maxPenalty <- sqrt(.Machine$double.xmax)

  ## For nearest neighbour bandwidths override default bandwidth.min
  ## and bandwidth.max and use sample size information.

  num.bw <- ncol(xdat)
  num.obs <- nrow(xdat)
  
  if(bwtype!="fixed") {
    bandwidth.min <- 1
    bandwidth.max <- num.obs-1
  }

  xdat <- as.data.frame(xdat)
  
  if(!is.null(bandwidth) && (length(bandwidth) != num.bw)) stop(" Error: bandwidth supplied but length not compatible with X data")

  if(is.null(nmulti)) nmulti <- min(5,num.bw)

  ## Determine which predictors are categorical and which are
  ## discrete... we care about unordered categorical kernels if the
  ## Aitchison & Aitken kernel is used since its bandwidth bounds are
  ## [0,(c-1)/c] and not 0/1 as are the rest of the unordered and
  ## ordered kernel bandwidth bounds.

  xdat.numeric <- sapply(1:num.bw,function(i){is.numeric(xdat[,i])})
  num.numeric <- ncol(as.data.frame(xdat[,xdat.numeric]))

  xdat.unordered <- sapply(1:num.bw,function(i){is.factor(xdat[,i])&&!is.ordered(xdat[,i])})
  num.unordered <- ncol(as.data.frame(xdat[,xdat.unordered]))

  if(numeric.scale==TRUE) {
    xdat.numeric.orig <- xdat[,xdat.numeric]
    xdat.scale <- attr(scale(xdat[,xdat.numeric]),"scaled:scale")
    xdat[,xdat.numeric] <- scale(xdat[,xdat.numeric])
    ## Note that scale() returns a matrix as the variable when there
    ## is only one variable. Check for this case and set variable to
    ## type numeric when this occurs.
    if(num.numeric==1) xdat[,xdat.numeric] <- as.numeric(xdat[,xdat.numeric])
  }

  ## Use degree for initial values if provided

  if(is.null(degree)) {
    if(cv == "degree-bandwidth") {
      degree <- sample(degree.min:degree.max, num.numeric, replace=T)
    }
    else {
      stop(paste(" Error: degree must be given when optimizing only bandwidth"))
    }
  }

  ## Use bandwidth for initial values if provided

  if(is.null(bandwidth)) {
    init.search.vals <- runif(num.bw,0,1)
    for(i in 1:num.bw) {
      if(xdat.numeric[i]==TRUE && bwtype=="fixed") {
        init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
      } 
      if(xdat.numeric[i]==TRUE && bwtype!="fixed") {
        init.search.vals[i] <- round(runif(1,2,sqrt(num.obs)))
      }
      if(xdat.unordered[i]==TRUE && ukertype=="aitchisonaitken") {
        c.num <- length(unique(xdat[,i]))
        init.search.vals[i] <- runif(1,0,(c.num-1)/c.num)
      }
    }
  } else {
    init.search.vals <- bandwidth
  }

  ## Create the function wrappers to be fed to the snomadr solver for
  ## leave-one-out cross-validation and Hurvich, Simonoff, and Tsai's
  ## AIC_c approach

  eval.lscv <- function(input, params){

    ydat <- params$ydat
    xdat <- params$xdat
    xdat.numeric <- params$xdat.numeric
    num.bw <- params$num.bw
    num.numeric <- params$num.numeric
    maxPenalty <- params$maxPenalty
    degree <- params$degree
    cv <- params$cv
    ukertype <- params$ukertype
    okertype <- params$okertype
    bwtype <- params$bwtype    

    bw.gamma <- input[1:num.bw]
    if(cv=="degree-bandwidth")
      degree <- round(input[(num.bw+1):(num.bw+num.numeric)])

    W <- W.glp(xdat,degree)

    if(all(bw.gamma>=0)&&all(bw.gamma[!xdat.numeric]<=1)) {
      lscv <- minimand.cv.ls(bws=bw.gamma,
                             ydat=ydat,
                             xdat=xdat,
                             degree=degree,
                             W=W,
                             ukertype=ukertype,
                             okertype=okertype,
                             bwtype=bwtype,                    
                             ...)
    } else {
      lscv <- maxPenalty
    }

    return(lscv)
  }

  eval.aicc <- function(input, params){

    ydat <- params$ydat
    xdat <- params$xdat
    xdat.numeric <- params$xdat.numeric
    num.bw <- params$num.bw
    num.numeric <- params$num.numeric
    maxPenalty <- params$maxPenalty
    degree <- params$degree
    cv <- params$cv
    ukertype <- params$ukertype
    okertype <- params$okertype
    bwtype <- params$bwtype

    bw.gamma <- input[1:num.bw]
    if(cv=="degree-bandwidth")
      degree <- round(input[(num.bw+1):(num.bw+num.numeric)])

    W <- W.glp(xdat,degree)

    if(all(bw.gamma>=0)&&all(bw.gamma[!xdat.numeric]<=1)) {
      aicc <- minimand.cv.aic(bws=bw.gamma,
                              ydat=ydat,
                              xdat=xdat,
                              degree=degree,
                              W=W,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,                    
                              ...)
    } else {
      aicc <- maxPenalty
    }

    return(aicc)
  }

  ## Generate the params fed to the snomadr solver

  params <- list()
  params$xdat <- xdat
  params$ydat <- ydat
  params$xdat.numeric <- xdat.numeric
  params$num.bw <- num.bw
  params$num.numeric <- num.numeric
  params$maxPenalty <- maxPenalty
  params$cv <- cv
  params$degree <- degree
  params$ukertype <- ukertype
  params$okertype <- okertype
  params$bwtype <- bwtype

  if(cv=="degree-bandwidth") {
    bbin <- c(rep(0, num.bw), rep(1, num.numeric))
    lb <- c(rep(bandwidth.min, num.bw), rep(degree.min, num.numeric))
    ub <- c(rep(bandwidth.max, num.bw), rep(degree.max, num.numeric))
  } else {
    bbin <- c(rep(0, num.bw))
    lb <- c(rep(bandwidth.min, num.bw))
    ub <- c(rep(bandwidth.max, num.bw))
  }

  for(i in 1:num.bw) {
    ## Need to do integer search for numeric predictors when bwtype is
    ## a nearest-neighbour, so set bbin appropriately.
    if(xdat.numeric[i]==TRUE && bwtype!="fixed") {
      bbin[i] <- 1
    }
    if(xdat.numeric[i]!=TRUE) {
      lb[i] <- 0.0
      ub[i] <- 1.0
    }
    ## Check for unordered and Aitchison/Aitken kernel
    if(xdat.unordered[i]==TRUE && ukertype=="aitchisonaitken") {
      c.num <- length(unique(xdat[,i]))
      ub[i] <- (c.num-1)/c.num
    }
  }

  ## No constraints

  bbout <-c(0)

  ## Multistarting

  fv.vec <- numeric(nmulti)

  ## Whether or not to display the information in snomadr
  print.output <- FALSE
  console <- newLineConsole()
  if(!is.null(opts$DISPLAY_DEGREE)){
    if(opts$DISPLAY_DEGREE>0){
      print.output <-TRUE
      console <- printPush("\rCalling NOMAD (Nonsmooth Optimization by Mesh Adaptive Direct Search)\n",console = console)
    }
  } else {
    print.output <-TRUE
    console <- printPush("\rCalling NOMAD (Nonsmooth Optimization by Mesh Adaptive Direct Search)\n",console = console)
  }

  degree.opt <- degree

  ## Generate all initial points for the multiple restarting
	x0.pts <- matrix(numeric(1), nmulti, length(bbin))
	for(iMulti in 1:nmulti) {
    ## First initialize to values for factors (`liracine' kernel)
    if(iMulti != 1) {
      init.search.vals <- runif(num.bw,0,1)
      for(i in 1:num.bw) {
        if(xdat.numeric[i]==TRUE && bwtype=="fixed") {
          init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
        }
        if(xdat.numeric[i]==TRUE && bwtype!="fixed") {
          init.search.vals[i] <- round(runif(1,2,sqrt(num.obs)))
        }
        if(xdat.unordered[i]==TRUE && ukertype=="aitchisonaitken") {
          c.num <- length(unique(xdat[,i]))
          init.search.vals[i] <- runif(1,0,(c.num-1)/c.num)
        }
      }
    }
    
    if(cv == "degree-bandwidth" && iMulti != 1)
      degree <- sample(degree.min:degree.max, num.numeric, replace=T)
    
    if(cv =="degree-bandwidth"){
      x0.pts[iMulti, ] <- c(init.search.vals, degree)
		}
    else {
      x0.pts[iMulti, ] <- c(init.search.vals)
		}
    
	}
	if(bwmethod == "cv.ls" ) {
			solution<-snomadr(eval.f=eval.lscv,
												n=length(bbin),
												x0=as.numeric(x0.pts),
												bbin=bbin,
												bbout=bbout,
												lb=lb,
												ub=ub,
												nmulti=nmulti,
												random.seed=random.seed,
												opts=opts,
												print.output=print.output,
												params=params);

	} else {
			solution<-snomadr(eval.f=eval.aicc,
												n=length(bbin),
												x0=as.numeric(x0.pts),
												bbin=bbin,
												bbout=bbout,
												lb=lb,
												ub=ub,
												nmulti=nmulti,
												random.seed=random.seed,
												opts=opts,
												print.output=print.output,
												params=params);
	}

	fv.vec[1] <- solution$objective

	bw.opt <- solution$solution[1:num.bw]
	if(numeric.scale==TRUE && bwtype=="fixed") {
			bw.opt[xdat.numeric] <- bw.opt[xdat.numeric]*xdat.scale
	}

	if(cv == "degree-bandwidth") {
			degree.opt <- solution$solution[(num.bw+1):(num.bw+num.numeric)]
	}
	fv <- solution$objective

	best <- NULL
	numimp <- 0   

  if(any(degree.opt==degree.max)) warning(paste(" an optimal degree equals search maximum (", degree.max,"): rerun with larger degree.max",sep=""))

  if(any(bw.opt==bandwidth.min)) warning(paste(" an optimal bandwidth equals search maximum (", bandwidth.min,"): rerun with smaller bandwidth.min",sep=""))

  console <- printPush("\r                        ",console = console)
  console <- printClear(console)
  console <- printPop(console)

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  return(list(bws=bw.opt,
              fv=fv,
              numimp=numimp,
              best=best,
              fv.vec=fv.vec,
              degree=degree.opt,
              bwtype=bwtype,
              ukertype=ukertype,
              okertype=okertype))

}

compute.bootstrap.errors <- function(tydat,
                                     txdat,
                                     exdat,
                                     eydat,
                                     bws,
                                     degree,
                                     ukertype,
                                     okertype,
                                     bwtype,boot.object=c("fitted","gradient","gradient.categorical"),
                                     plot.errors.boot.num=99,
                                     plot.errors.type=c("quantiles","standard"),
                                     plot.errors.quantiles=c(.025,.975),
                                     alpha=0.05,
                                     gradient.vec=NULL,
                                     gradient.categorical=FALSE,
                                     gradient.categorical.index=NULL,                                     
                                     ...){
  
  plot.errors.type <- match.arg(plot.errors.type)
  boot.object <- match.arg(boot.object)

  if(missing(exdat)) {
    neval <- nrow(txdat)
    exdat <- txdat
  } else {
    neval <- nrow(exdat)
  }

  boot.err <- matrix(data = NA, nrow = neval, ncol = 2)

  ## Conduct simple iid residual bootstrap

  boot.func.mean <- function(model.fitted,indices){
    est.boot <- npglpreg(tydat=model.fitted+(tydat-model.fitted)[indices],
                         txdat=txdat,
                         exdat=exdat,
                         bws=bws,
                         degree=degree,
                         ukertype=ukertype,
                         okertype=okertype,
                         bwtype=bwtype,
                         gradient.vec=gradient.vec,
                         gradient.categorical=gradient.categorical,
                         gradient.categorical.index=gradient.categorical.index,
                         ...)
    if(boot.object=="fitted") {
      return(est.boot$fitted.values)
    } else if(boot.object=="gradient") {
      return(est.boot$gradient)
    } else if(boot.object=="gradient.categorical") {
      return(est.boot$gradient.categorical.mat[,gradient.categorical.index])
    }
  }
  
  ## Fitted values for the training data required
  
  est <- npglpreg(tydat=tydat,
                  txdat=txdat,
                  bws=bws,
                  degree=degree,
                  ukertype=ukertype,
                  okertype=okertype,
                  bwtype=bwtype,
                  ...)
  
  model.fitted <- est$fitted.values
  
  boot.out <- boot(data = model.fitted,
                   statistic = boot.func.mean,
                   R = plot.errors.boot.num)
  
  if (plot.errors.type == "standard") {
    boot.err[,1:2] <- qnorm(1-alpha/2)*sqrt(diag(cov(boot.out$t)))
    boot.err[,1] <- boot.out$t0 - boot.err[,1]
    boot.err[,2] <- boot.out$t0 + boot.err[,2]
  }
  else if (plot.errors.type == "quantiles") {
    boot.err[,1:2] <- t(sapply(as.data.frame(boot.out$t),
                               function (y) {
                                 quantile(y,probs = plot.errors.quantiles)
                               }))
  }
  
  return(cbind(boot.out$t0,boot.err))

}

plot.npglpreg <- function(x,
                          mean=TRUE,
                          deriv=0,
                          ci=FALSE,
                          num.eval=100,
                          common.scale=TRUE,
                          xtrim = 0.0,
                          xq = 0.5,
                          plot.behavior = c("plot","plot-data","data"),
                          plot.errors.boot.num=99,
                          plot.errors.type=c("quantiles","standard"),
                          plot.errors.quantiles=c(.025,.975),
                          persp.rgl=FALSE,
                          ...) {
  
  plot.behavior <- match.arg(plot.behavior)
  plot.errors.type <- match.arg(plot.errors.type)

  object <- x

  ## Needed for correctly obtaining predictions

  degree <- object$degree
  bws <- object$bws
  bwtype <- object$bwtype
  ukertype <- object$ukertype
  okertype <- object$okertype
  
  txdat <- object$x
  tydat <- object$y

  xq <- double(ncol(txdat)) + xq

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("\rWorking...",console = console)

  ## Mean
  
  if(mean==TRUE && deriv==0) {

    if(!persp.rgl) {
      
      mg <- list()
      
      for(i in 1:NCOL(object$x)) {
        
        if(!is.factor(object$x[,i])) {
          exdat <- matrix(NA,nrow=num.eval,ncol=NCOL(object$x))
          neval <- num.eval
        } else {
          neval <- length(unique(object$x[,i]))
          exdat <- matrix(NA,nrow=neval,ncol=NCOL(object$x))
        }
        
        exdat <- data.frame(exdat)
        
        if(!is.factor(object$x[,i])) {
          xlim <- trim.quantiles(object$x[,i],xtrim)          
          exdat[,i] <- seq(xlim[1],xlim[2],length=neval)
        } else {
          exdat[,i] <- sort(unique(object$x[,i]))
        }
        
        for(j in (1:NCOL(object$x))[-i]) {
          exdat[,j] <- rep(uocquantile(object$x[,j],prob=xq[j]),neval)
        }

        names(exdat) <- object$xnames

        est <- npglpreg.default(tydat=tydat,
                                txdat=txdat,
                                exdat=exdat,
                                bws=bws,
                                degree=degree,
                                ukertype=ukertype,
                                okertype=okertype,
                                bwtype=bwtype,
                                ...)

        fitted.values <- est$fitted.values

        if(!ci) {
          
          mg[[i]] <- data.frame(exdat[,i],fitted.values)
          names(mg[[i]]) <- c(names(exdat)[i],"mean")
          
        } else {

          console <- printClear(console)
          console <- printPop(console)
          console <- printPush(paste("\rConducting ",plot.errors.boot.num," bootstrap resamples for predictor ",i,"...",sep=""),console = console)
          
          ci.out <- compute.bootstrap.errors(tydat=tydat,
                                             txdat=txdat,
                                             exdat=exdat,
                                             bws=bws,
                                             degree=degree,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,
                                             boot.object="fitted",
                                             plot.errors.boot.num=plot.errors.boot.num,
                                             plot.errors.type=plot.errors.type,
                                             plot.errors.quantiles=plot.errors.quantiles)
          
          mg[[i]] <- data.frame(exdat[,i],ci.out)
          names(mg[[i]]) <- c(names(exdat)[i],"mean","lwr","upr")
          
        }

      }
      
      if(common.scale) {
        min.mg <- Inf
        max.mg <- -Inf
        for(i in 1:length(mg)) {
          min.mg <- min(min.mg,min(mg[[i]][,-1]))
          max.mg <- max(max.mg,max(mg[[i]][,-1]))
        }
        ylim <- c(min.mg,max.mg)
      } else {
        ylim <- NULL
      }
      
      if(plot.behavior!="data") {

        if(!is.null(object$num.categorical)||(object$num.numeric>1)) par(mfrow=dim.plot(NCOL(object$x)))
      
        for(i in 1:NCOL(object$x)) {

          if(!ci) {
        
            plot(mg[[i]][,1],mg[[i]][,2],
                 xlab=names(exdat)[i],
                 ylab="Conditional Mean",
                 ylim=ylim,
                 type="l")
            
          } else {
            if(!common.scale) ylim <- c(min(mg[[i]][,-1]),max(mg[[i]][,-1]))
            plot(mg[[i]][,1],mg[[i]][,2],
                 xlab=names(exdat)[i],
                 ylab="Conditional Mean",
                 ylim=ylim,
                 type="l")
            ## Need to overlay for proper plotting of factor errorbars
            par(new=TRUE)
            plot(mg[[i]][,1],mg[[i]][,3],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2)
            par(new=TRUE)
            plot(mg[[i]][,1],mg[[i]][,4],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2)
          }
          
        }
        
      }
      
    } else {
      
      if(!require(rgl)) stop(" Error: you must first install the rgl package")
      
      if(object$num.categorical != 0) stop(" Error: persp3d is for continuous predictors only")
      if(object$num.numeric != 2) stop(" Error: persp3d is for cases involving two continuous predictors only")
      
      newdata <- matrix(NA,nrow=num.eval,ncol=2)
      newdata <- data.frame(newdata)
      
      xlim <- trim.quantiles(object$x[,1],xtrim)
      ylim <- trim.quantiles(object$x[,2],xtrim)      
      
      x1.seq <- seq(xlim[1],xlim[2],length=num.eval)
      x2.seq <- seq(ylim[1],ylim[2],length=num.eval)    
      x.grid <- expand.grid(x1.seq,x2.seq)
      newdata <- data.frame(x.grid[,1],x.grid[,2])
      names(newdata) <- names(object$x)
      
      z <- matrix(predict(object,newdata=newdata),num.eval,num.eval)

      mg <- list()

      mg[[1]] <- data.frame(newdata,z)

      if(plot.behavior!="data") {
        
        num.colors <- 1000
        colorlut <- topo.colors(num.colors) 
        col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]
        
        open3d()
        
        par3d(windowRect=c(900,100,900+640,100+640))
        rgl.viewpoint(theta = 0, phi = -70, fov = 80)
        
        persp3d(x=x1.seq,y=x2.seq,z=z,
                xlab=names(object$x)[1],ylab=names(object$x)[2],zlab="Y",
                ticktype="detailed",      
                border="red",
                color=col,
                alpha=.7,
                back="lines",
                main="Conditional Mean")
        
        grid3d(c("x", "y+", "z"))
        
        play3d(spin3d(axis=c(0,0,1), rpm=5), duration=15)
        
      }

    }
    
    if(plot.behavior!="plot") {
      console <- printClear(console)
      console <- printPop(console)
      return(mg)
    }
      
  }
    
  ## deriv

  if(deriv > 0) {
    
    rg <- list()

    i.numeric <- 1
    i.categorical <- 1
    
    for(i in 1:NCOL(object$x)) {

      gradient.vec <- NULL
      
      if(!is.factor(object$x[,i])) {
        newdata <- matrix(NA,nrow=num.eval,ncol=NCOL(object$x))
        neval <- num.eval
        gradient.vec <- rep(0,object$num.numeric)
        gradient.vec[i.numeric] <- deriv
      } else {
        neval <- length(unique(object$x[,i]))
        newdata <- matrix(NA,nrow=neval,ncol=NCOL(object$x))
      }
      
      newdata <- data.frame(newdata)
      
      if(!is.factor(object$x[,i])) {
        xlim <- trim.quantiles(object$x[,i],xtrim)        
        newdata[,i] <- seq(xlim[1],xlim[2],length=neval)
      } else {
        newdata[,i] <- sort(unique(object$x[,i]))
      }
      
      for(j in (1:NCOL(object$x))[-i]) {
        newdata[,j] <- rep(uocquantile(object$x[,j],prob=xq[j]),neval)
      }
      
      names(newdata) <- object$xnames

      est <- npglpreg.default(tydat=tydat,
                              txdat=txdat,
                              exdat=newdata,
                              bws=bws,
                              degree=degree,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,
                              gradient.vec=gradient.vec,
                              gradient.categorical=TRUE,
                              ...)

      if(!is.factor(object$x[,i])) {
        fitted.values <- est$gradient
      } else {
        fitted.values <- est$gradient.categorical.mat[,i.categorical]
      }
        

      if(!ci) {
        
        rg[[i]] <- data.frame(newdata[,i],fitted.values)
        names(rg[[i]]) <- c(names(newdata)[i],"deriv")
        
      } else {
        
        console <- printClear(console)
        console <- printPop(console)
        console <- printPush(paste("\rConducting ",plot.errors.boot.num," bootstrap resamples for predictor ",i,"...",sep=""),console = console)
          
        if(!is.factor(object$x[,i])) {
          ci.out <- compute.bootstrap.errors(tydat=tydat,
                                             txdat=txdat,
                                             exdat=newdata,
                                             bws=bws,
                                             degree=degree,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,
                                             boot.object="gradient",
                                             plot.errors.boot.num=plot.errors.boot.num,
                                             plot.errors.type=plot.errors.type,
                                             plot.errors.quantiles=plot.errors.quantiles,
                                             gradient.vec=gradient.vec)
        } else {
          ci.out <- compute.bootstrap.errors(tydat=tydat,
                                             txdat=txdat,
                                             exdat=newdata,
                                             bws=bws,
                                             degree=degree,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,
                                             boot.object="gradient.categorical",
                                             plot.errors.boot.num=plot.errors.boot.num,
                                             plot.errors.type=plot.errors.type,
                                             plot.errors.quantiles=plot.errors.quantiles,
                                             gradient.categorical=TRUE,                                
                                             gradient.categorical.index=i.categorical)
        }
        
        rg[[i]] <- data.frame(newdata[,i],ci.out)
        names(rg[[i]]) <- c(names(newdata)[i],"deriv","lwr","upr")
        
      }
      
      if(!is.factor(object$x[,i])) {
        i.numeric <- i.numeric + 1
      } else {
        i.categorical <- i.categorical + 1
      }          
        
    }
    
    if(common.scale) {
      min.rg <- Inf
      max.rg <- -Inf
      for(i in 1:length(rg)) {
        min.rg <- min(min.rg,min(rg[[i]][,-1]))
        max.rg <- max(max.rg,max(rg[[i]][,-1]))
      }
      ylim <- c(min.rg,max.rg)
    } else {
      ylim <- NULL
    }
    
    if(plot.behavior!="data") {
      
      if(!is.null(object$num.categorical)||(object$num.numeric>1)) par(mfrow=dim.plot(NCOL(object$x)))
    
      for(i in 1:NCOL(object$x)) {
        
          if(!ci) {
          plot(rg[[i]][,1],rg[[i]][,2],
               xlab=names(newdata)[i],
               ylab=ifelse(!is.factor(newdata[,i]), paste("Order", deriv,"Gradient"), "Difference in Levels"),
               ylim=ylim,
               type="l")
          
        } else {
          if(!common.scale) ylim <- c(min(rg[[i]][,-1]),max(rg[[i]][,-1]))
          plot(rg[[i]][,1],rg[[i]][,2],
               xlab=names(newdata)[i],
               ylab=ifelse(!is.factor(newdata[,i]), paste("Order", deriv,"Gradient"), "Difference in Levels"),
               ylim=ylim,
               type="l")
          ## Need to overlay for proper plotting of factor errorbars
          par(new=TRUE)
          plot(rg[[i]][,1],rg[[i]][,3],
               xlab="",
               ylab="",
               ylim=ylim,
               type="l",
               axes=FALSE,
               lty=2,
               col=2)
          par(new=TRUE)
          plot(rg[[i]][,1],rg[[i]][,4],
               xlab="",
               ylab="",
               ylim=ylim,
               type="l",
               axes=FALSE,
               lty=2,
               col=2)
        }
        
      }
      
    }
    
    if(plot.behavior!="plot") {
      console <- printClear(console)
      console <- printPop(console)
      return(rg)
    }
    
  }

  console <- printClear(console)
  console <- printPop(console)

  ## Reset par to 1,1 (can be modified above)
  
  if(!persp.rgl) par(mfrow=c(1,1))
  
}
