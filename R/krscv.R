## This function conducts kernel regression spline
## cross-validation. It takes as input a data.frame xz containing a
## mix of numeric and factor predictors and a vector y. A range of
## arguments can be provided, and one can do search on the bandwidths
## and both the degree and knots ("degree-knots") or the degree
## holding the number of knots (segments+1) constant or the number of
## knots (segments+1) holding the degree constant. A variety of basis
## types are supported (functional anova "additive-tensor",
## "additive", or "tensor") and the argument "auto" will evaluate
## choose the basis type automatically.

## Currently search is exhaustive taking basis.maxdim as the maximum
## number of the spline degree (0,1,...) and number of segments
## (1,2,...). This is a quadratic integer programming problem so
## ideally I require an IQP (MIQP for kernel-weighting)
## solver. Currently in R there is no such beast.

krscv <- function(xz,
                  y,
                  basis.maxdim=5,
                  kernel.type=c("nominal","ordinal"),
                  restarts=0,
                  complexity=c("degree-knots","degree","knots"),
                  knots=c("quantiles","uniform"),
                  basis = c("additive-tensor","additive","tensor","auto"),
                  cv.norm=c("L2","L1"),
                  degree=degree,
                  segments=segments) {

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)
  cv.norm <- match.arg(cv.norm)  

  ## First define the cv function to be fed to optim

  kernel.type <- match.arg(kernel.type)

  t1 <- Sys.time()

  cv.func <- function(input,
                      x,
                      y,
                      z,
                      K,
                      basis.maxdim,
                      restart,
                      num.restarts,
                      z.unique,
                      ind,
                      ind.vals,
                      nrow.z.unique,
                      kernel.type,
                      j=NULL,
                      nrow.K.mat=NULL,
                      t2=NULL,
                      complexity=complexity,
                      knots=knots,
                      basis=basis,
                      cv.norm=cv.norm) {

    ## K is a matrix, column 1 degree, column 2 segments, either or
    ## both can be determined via cv so need to take care to allow
    ## user to select knots (degree fixed), degree (knots fixed), or
    ## both degree and knots. The values used to evaluate the cv
    ## function are passed below.
    
    if(is.null(K)) {
      num.x <- NCOL(x)
      num.z <- NCOL(z)
      K <- round(cbind(input[1:num.x],input[(num.x+1):(2*num.x)]))
      lambda <- input[(2*num.x+1):(2*num.x+num.z)]
    } else {
      K <- round(cbind(K[1:num.x],K[(num.x+1):(2*num.x)]))      
      lambda <- input
    }
    ## When using weights= lambda of zero fails. Trivial to trap.
    lambda <- ifelse(lambda <= 0, .Machine$double.eps, lambda)

    cv <- cv.kernel.spline(x=x,
                           y=y,
                           z=z,
                           K=K,
                           lambda=lambda,
                           z.unique=z.unique,
                           ind=ind,
                           ind.vals=ind.vals,
                           nrow.z.unique=nrow.z.unique,
                           kernel.type=kernel.type,
                           knots=knots,
                           basis=basis)

    ## Some i/o unless options(crs.messages=FALSE)

    fw.format.3 <- function(input) sapply(input,sprintf,fmt="%#.3f")
    fw.format.2 <- function(input) sapply(input,sprintf,fmt="%#.2f")

    if(complexity=="degree") {
      if(!is.null(j)) {
        if(j==1) {
          tmp.1 <- paste(j,"/",nrow.K.mat,", d[1]=",K[1,1],sep="")
        } else {
          dt <- (t2-t1)*(nrow.K.mat-j+1)/j
          tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                         fw.format.2(as.numeric((t2-t1),units="mins")),
                         "m",sep="")
          tmp.1 <- paste(j,"/",nrow.K.mat,tmp.0,", d[1]=",K[1,1],sep="")
        }
      } else {
        tmp.1 <- paste("d[1]=", K[1,1],sep="")
      }
      if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", d[", i, "]=", K[i,1],sep="")
    } else  if(complexity=="knots") {
      if(!is.null(j)) {
        if(j==1) {
          tmp.1 <- paste(j,"/",nrow.K.mat,", s[1]=",K[1,2],sep="")
        } else {
          dt <- (t2-t1)*(nrow.K.mat-j+1)/j
          tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                         fw.format.2(as.numeric((t2-t1),units="mins")),
                         "m",sep="")
          tmp.1 <- paste(j,"/",nrow.K.mat,tmp.0,", s[1]=",K[1,2],sep="")
        }
      } else {
        tmp.1 <- paste("s[1]=", K[1,2],sep="")
      }
      if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", s[", i, "]=", K[i,2],sep="")
    } else  if(complexity=="degree-knots") {
      if(!is.null(j)) {
        if(j==1) {
          tmp.1 <- paste(j,"/",nrow.K.mat,", d[1]=",K[1,1],sep="")
        } else {
          dt <- (t2-t1)*(nrow.K.mat-j+1)/j
          tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                         fw.format.2(as.numeric((t2-t1),units="mins")),
                         "m",sep="")
          tmp.1 <- paste(j,"/",nrow.K.mat,tmp.0,", d[1]=",K[1,1],sep="")
        }
      } else {
        tmp.1 <- paste("d[1]=", K[1,1],sep="")
      }
      if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", d[", i, "]=", K[i,1],sep="")
      for(i in 1:num.x) tmp.1 <- paste(tmp.1, ", s[", i, "]=", K[i,2],sep="")      
    }

    ## For i/o for z variables
    
    tmp.2 <- paste(", rs=", restart, "/", num.restarts,sep="")
    tmp.3 <- ""
    for(i in 1:num.z) tmp.3 <- paste(tmp.3, ", l[", i, "]=", fw.format.3(lambda[i]),sep="")
    tmp.4 <- paste(", cv=", format(cv,digits=6), sep="")
    if(num.restarts > 0) {
      msg <- paste(tmp.1,tmp.2,tmp.3,tmp.4,sep="")
    } else {
      msg <- paste(tmp.1,tmp.3,tmp.4,sep="")
    }

    console <<- printClear(console)
    console <<- printPush(msg,console = console)

    return(cv)

  }

  xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
  x <- xztmp$x
  z <- xztmp$z
  if(is.null(z)) stop(" categorical kernel smoothing requires ordinal/nominal predictors")

  z <- as.matrix(xztmp$z)
  num.z <- NCOL(z)
  z.unique <- uniquecombs(z)
  ind <-  attr(z.unique,"index")
  ind.vals <-  unique(ind) #sort(unique(ind))
  nrow.z.unique <- NROW(z.unique)
  num.x <- NCOL(x)
  n <- NROW(x)

  if(complexity=="degree") {
    if(missing(segments)) stop("segments missing for cross-validation of spline degree")
    if(length(segments)!=num.x) stop(" segments vector must be the same length as x")  
  } else if(complexity=="knots") {
    if(missing(degree)) stop("degree missing for cross-validation of number of spline knots")
    if(length(degree)!=num.x) stop(" degree vector must be the same length as x")
  }

  ## For kernel regression spline, if there is only one continuous
  ## predictor (i.e. num.x==1) disable auto, set to additive (which is
  ## additive-tensor and tensor in this case, so don't waste time
  ## doing all three).

  if(num.x==1 & basis == "auto") basis <- "additive"

  ## This call will trap zero/negative degrees of freedom immediately
  ## rather than letting cv proceed only to be halted after the lower
  ## dimension models have been estimated when this occurs.

##  if(complexity=="degree") {
##    if(basis == "auto") {
##      k <- max(c(ncol(prod.spline(x=x,K=cbind(rep(basis.maxdim,num.x),segments),knots=knots,basis="additive-tensor")),
##                 ncol(prod.spline(x=x,K=cbind(rep(basis.maxdim,num.x),segments),knots=knots,basis="additive")),
##                 ncol(prod.spline(x=x,K=cbind(rep(basis.maxdim,num.x),segments),knots=knots,basis="tensor"))))
##    } else {
##      k <- ncol(prod.spline(x=x,K=cbind(rep(basis.maxdim,num.x),segments),knots=knots,basis=basis))
##    }
##  } else if(complexity=="knots") {
##    if(basis == "auto") {
##      k <- max(c(ncol(prod.spline(x=x,K=cbind(degree,rep(basis.maxdim,num.x)),knots=knots,basis="additive-tensor")),
##                 ncol(prod.spline(x=x,K=cbind(degree,rep(basis.maxdim,num.x)),knots=knots,basis="additive")),
##                 ncol(prod.spline(x=x,K=cbind(degree,rep(basis.maxdim,num.x)),knots=knots,basis="tensor"))))
##    } else {
##      k <- ncol(prod.spline(x=x,K=cbind(degree,rep(basis.maxdim,num.x)),knots=knots,basis=basis))
##    }
##  } else if(complexity=="degree-knots"){
##    if(basis == "auto") {
##      k <- max(c(ncol(prod.spline(x=x,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),knots=knots,basis="additive-tensor")),
##                 ncol(prod.spline(x=x,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),knots=knots,basis="additive")),
##                 ncol(prod.spline(x=x,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),knots=knots,basis="tensor"))))
##    } else {
##      k <- ncol(prod.spline(x=x,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),knots=knots,basis=basis))
##    }
##  }
##
##  df <- n - k
##
##  if(df <= 0) {
##    stop(paste(" maximum basis dimension (",k,") would equal/exceed sample size (",n,")\n   perhaps use basis=\"additive\" or else decrease basis.maxdim",sep=""))
##  } else if(df <= 10) {
##    warning(paste(" maximum basis dimension (",k,") and sample size (",n,") close",sep=""))
##  }

# 17/06/10 - really dig into what this tests and whether it is necessary or not (probably is!)
#  if(min(table(ind)) <= k) stop(paste(" insufficient data for one or more unique combinations of z (",min(table(ind))," obs.)\n   in order to estimate spline at basis.maxdim (",k," bases):\n   either reduce basis.maxdim or collapse categories",sep=""))

  if(basis.maxdim < 1) stop(" basis.maxdim must be greater than or equal to 1")

  console <- newLineConsole()
  console <- printPush("Working...",console = console)

  ## Exhaustive evaluation over all combinations of K, search over
  ## lambda for each combination

  if(complexity!="degree-knots") {
    K.mat <- matrix.combn(0:basis.maxdim,num.x)
  } else {
    K.mat <- matrix.combn(0:basis.maxdim,2*num.x)
  }

  nrow.K.mat <- NROW(K.mat)
  cv.vec <- numeric(nrow.K.mat)
  basis.vec <- character(nrow.K.mat)
  lambda.mat <- matrix(NA,nrow.K.mat,num.z)

  t2 <- Sys.time() ## placeholder

  output <- list()
  output.restart <- list()

  ## First need to create KI mat outside loop and not - back to cbind issue...
  
  if(complexity=="degree") {
    K.mat <- cbind(K.mat[,1:num.x],matrix(segments,nrow(K.mat),length(segments),byrow=TRUE))
  } else if(complexity=="knots") {
    K.mat <- cbind(matrix(degree,nrow(K.mat),length(degree),byrow=TRUE),K.mat[,1:num.x]+1) 
  } else if(complexity=="degree-knots") {
    K.mat[,(num.x+1):(2*num.x)] <- K.mat[,(num.x+1):(2*num.x)]+1
  }

  for(j in 1:nrow.K.mat) {

    ## Initialize    

    cv.vec[j] <- .Machine$double.xmax

    if(basis=="auto") {

      ## First basis=="additive-tensor"

      output$convergence <- 42

      while(output$convergence != 0) {
        output <- optim(par=runif(num.z),
                        cv.func,
                        lower=rep(0,num.z),
                        upper=rep(1,num.z),
                        method="L-BFGS-B",
                        x=x,
                        y=y,
                        z=z,
                        K=K.mat[j,],
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=restarts,
                        z.unique=z.unique,
                        ind=ind,
                        ind.vals=ind.vals,
                        nrow.z.unique=nrow.z.unique,
                        kernel.type=kernel.type,
                        j=j,
                        nrow.K.mat=nrow.K.mat,
                        t2=t2,
                        complexity=complexity,
                        knots=knots,
                        basis="additive-tensor")

      }

      if(restarts > 0) {

        for(r in 1:restarts) {
          
          output.restart$convergence <- 42
          
          while(output.restart$convergence != 0) {

            output.restart <- optim(par=runif(num.z),
                                    cv.func,
                                    lower=rep(0,num.z),
                                    upper=rep(1,num.z),
                                    method="L-BFGS-B",
                                    x=x,
                                    y=y,
                                    z=z,
                                    K=K.mat[j,],
                                    basis.maxdim=basis.maxdim,
                                    restart=r,
                                    num.restarts=restarts,
                                    z.unique=z.unique,
                                    ind=ind,
                                    ind.vals=ind.vals,
                                    nrow.z.unique=nrow.z.unique,
                                    kernel.type=kernel.type,
                                    j=j,
                                    nrow.K.mat=nrow.K.mat,
                                    t2=t2,
                                    complexity=complexity,
                                    knots=knots,
                                    basis="additive-tensor")

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      if(output$value < cv.vec[j]) {
        cv.vec[j] <- output$value
        basis.vec[j] <- "additive-tensor"
        lambda.mat[j,] <- output$par
      }

      ## Next, basis=="additive"

      output$convergence <- 42

      while(output$convergence != 0) {

        output <- optim(par=runif(num.z),
                        cv.func,
                        lower=rep(0,num.z),
                        upper=rep(1,num.z),
                        method="L-BFGS-B",
                        x=x,
                        y=y,
                        z=z,
                        K=K.mat[j,],
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=restarts,
                        z.unique=z.unique,
                        ind=ind,
                        ind.vals=ind.vals,
                        nrow.z.unique=nrow.z.unique,
                        kernel.type=kernel.type,
                        j=j,
                        nrow.K.mat=nrow.K.mat,
                        t2=t2,
                        complexity=complexity,
                        knots=knots,
                        basis="additive")

      }

      if(restarts > 0) {

        for(r in 1:restarts) {

          output.restart$convergence <- 42

          while(output.restart$convergence != 0) {

            output.restart <- optim(par=runif(num.z),
                                    cv.func,
                                    lower=rep(0,num.z),
                                    upper=rep(1,num.z),
                                    method="L-BFGS-B",
                                    x=x,
                                    y=y,
                                    z=z,
                                    K=K.mat[j,],
                                    basis.maxdim=basis.maxdim,
                                    restart=r,
                                    num.restarts=restarts,
                                    z.unique=z.unique,
                                    ind=ind,
                                    ind.vals=ind.vals,
                                    nrow.z.unique=nrow.z.unique,
                                    kernel.type=kernel.type,
                                    j=j,
                                    nrow.K.mat=nrow.K.mat,
                                    t2=t2,
                                    complexity=complexity,
                                    knots=knots,
                                    basis="additive")

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      if(output$value < cv.vec[j]) {
        cv.vec[j] <- output$value
        basis.vec[j] <- "additive"
        lambda.mat[j,] <- output$par
      }

      ## Next, basis=="tensor"

      output$convergence <- 42

      while(output$convergence != 0) {

        output <- optim(par=runif(num.z),
                        cv.func,
                        lower=rep(0,num.z),
                        upper=rep(1,num.z),
                        method="L-BFGS-B",
                        x=x,
                        y=y,
                        z=z,
                        K=K.mat[j,],
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=restarts,
                        z.unique=z.unique,
                        ind=ind,
                        ind.vals=ind.vals,
                        nrow.z.unique=nrow.z.unique,
                        kernel.type=kernel.type,
                        j=j,
                        nrow.K.mat=nrow.K.mat,
                        t2=t2,
                        complexity=complexity,
                        knots=knots,
                        basis="tensor")

      }

      if(restarts > 0) {

        for(r in 1:restarts) {

          output.restart$convergence <- 42

          while(output.restart$convergence != 0) {

            output.restart <- optim(par=runif(num.z),
                                    cv.func,
                                    lower=rep(0,num.z),
                                    upper=rep(1,num.z),
                                    method="L-BFGS-B",
                                    x=x,
                                    y=y,
                                    z=z,
                                    K=K.mat[j,],
                                    basis.maxdim=basis.maxdim,
                                    restart=r,
                                    num.restarts=restarts,
                                    z.unique=z.unique,
                                    ind=ind,
                                    ind.vals=ind.vals,
                                    nrow.z.unique=nrow.z.unique,
                                    kernel.type=kernel.type,
                                    j=j,
                                    nrow.K.mat=nrow.K.mat,
                                    t2=t2,
                                    complexity=complexity,
                                    knots=knots,
                                    basis="tensor")

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      if(output$value < cv.vec[j]) {
        cv.vec[j] <- output$value
        basis.vec[j] <- "tensor"
        lambda.mat[j,] <- output$par
      }

    } else { ## end auto

      ## Either basis=="additive-tensor" or "additive" or "tensor"

      output$convergence <- 42

      while(output$convergence != 0) {

        output <- optim(par=runif(num.z),
                        cv.func,
                        lower=rep(0,num.z),
                        upper=rep(1,num.z),
                        method="L-BFGS-B",
                        x=x,
                        y=y,
                        z=z,
                        K=K.mat[j,],
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=restarts,
                        z.unique=z.unique,
                        ind=ind,
                        ind.vals=ind.vals,
                        nrow.z.unique=nrow.z.unique,
                        kernel.type=kernel.type,
                        j=j,
                        nrow.K.mat=nrow.K.mat,
                        t2=t2,
                        complexity=complexity,
                        knots=knots,
                        basis=basis)

      }

      if(restarts > 0) {

        for(r in 1:restarts) {

          output.restart$convergence <- 42

          while(output.restart$convergence != 0) {

            output.restart <- optim(par=runif(num.z),
                                    cv.func,
                                    lower=rep(0,num.z),
                                    upper=rep(1,num.z),
                                    method="L-BFGS-B",
                                    x=x,
                                    y=y,
                                    z=z,
                                    K=K.mat[j,],
                                    basis.maxdim=basis.maxdim,
                                    restart=r,
                                    num.restarts=restarts,
                                    z.unique=z.unique,
                                    ind=ind,
                                    ind.vals=ind.vals,
                                    nrow.z.unique=nrow.z.unique,
                                    kernel.type=kernel.type,
                                    j=j,
                                    nrow.K.mat=nrow.K.mat,
                                    t2=t2,
                                    complexity=complexity,
                                    knots=knots,
                                    basis=basis)

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      if(output$value < cv.vec[j]) {
        cv.vec[j] <- output$value
        basis.vec[j] <- basis
        lambda.mat[j,] <- output$par
      }

    }

    t2 <- Sys.time()

  }

  ## Sort on cv.vec

  ocv.vec <- order(cv.vec)

  cv.min <- cv.vec[ocv.vec][1]
  if(complexity=="degree") {
    K.opt <- c(K.mat[ocv.vec,1:num.x,drop=FALSE][1,],segments)
  } else if(complexity=="knots") {
    K.opt <- c(degree,K.mat[ocv.vec,1:num.x,drop=FALSE][1,]+1) 
  } else if(complexity=="degree-knots") {
    K.opt <- K.mat[ocv.vec,,drop=FALSE][1,]
  }
  lambda.opt <- lambda.mat[ocv.vec,,drop=FALSE][1,]
  basis.opt <- basis.vec[ocv.vec][1]
  degree <- K.opt[1:num.x]
  segments <- K.opt[(num.x+1):(2*num.x)]
  if(!is.null(z)) I.opt <- K.opt[(2*num.x+1):(2*num.x+num.z)]
  
  console <- printClear(console)
  console <- printPop(console)

  if(any(degree==basis.maxdim)) warning(paste(" optimal degree equals search maximum (", basis.maxdim,"): rerun with larger basis.maxdim",sep=""))
  if(any(segments==(basis.maxdim+1))) warning(paste(" optimal segment equals search maximum (", basis.maxdim+1,"): rerun with larger basis.maxdim",sep=""))  

  crscv(K=K.opt,
        I=NULL,
        basis=basis.opt,
        basis.vec=basis.vec,
        basis.maxdim=basis.maxdim,
        complexity=complexity,
        knots=knots,
        degree=degree,
        segments=segments,
        restarts=restarts,
        K.mat=K.mat,
        lambda=lambda.opt,
        lambda.mat=lambda.mat,
        cv.func=cv.min,
        cv.func.vec=as.matrix(cv.vec))

}
