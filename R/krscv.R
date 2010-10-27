krscv <- function(xz,
                  y,
                  basis.maxdim=10,
                  kernel.type=c("nominal","ordinal"),
                  restarts=0,
                  complexity=c("degree","knots"),
                  basis = c("additive-tensor","additive","tensor","auto"),
                  cv.norm=c("L2","L1"),
                  degree=degree,
                  nbreak=nbreak) {

  complexity <- match.arg(complexity)
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
                      degree=degree,
                      nbreak=nbreak,
                      complexity=complexity,
                      basis=basis,
                      cv.norm=cv.norm) {

    ## For model of given complexity search for optimal bandwidths

    if(is.null(K)) {
      num.x <- NCOL(x)
      num.z <- NCOL(z)
      K <- round(input[1:num.x])
      lambda <- input[(num.x+1):(num.x+num.z)]
    } else {
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
                           degree=degree,
                           nbreak=nbreak,
                           complexity=complexity,
                           basis=basis)

    ## Some i/o unless options(crs.messages=FALSE)

    fw.format.3 <- function(input) sapply(input,sprintf,fmt="%#.3f")
    fw.format.2 <- function(input) sapply(input,sprintf,fmt="%#.2f")
    if(!is.null(j)) {
      if(j==1) {
        tmp.1 <- paste(j,"/",nrow.K.mat,", k[1]=",K[1],sep="")
      } else {
        dt <- (t2-t1)*(nrow.K.mat-j+1)/j
        tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                       fw.format.2(as.numeric((t2-t1),units="mins")),
                       "m",sep="")
        tmp.1 <- paste(j,"/",nrow.K.mat,tmp.0,", k[1]=",K[1],sep="")
      }
    } else {
      tmp.1 <- paste("k[1]=", K[1],sep="")
    }
    if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", k[", i, "]=", K[i],sep="")
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

  ## check whether tensor spline dimension results in zero/low df...

  k <- basis.maxdim^num.x + ifelse(basis!="additive",basis.maxdim*num.x,0)
  df <- n - k
  if(df <= 0) {
    stop(paste(" maximum spline dimension (",k,") equals/exceeds sample size (",n,")\n   perhaps use basis=\"additive\" or else decrease basis.maxdim",sep=""))
  } else if(df <= 10) {
    warning(paste(" maximum spline dimension (",k,") and sample size (",n,") close",sep=""))
  }
# 17/06/10 - really dig into what this tests and whether it is necessary or not (probably is!)
#  if(min(table(ind)) <= k) stop(paste(" insufficient data for one or more unique combinations of z (",min(table(ind))," obs.)\n   in order to estimate spline at basis.maxdim (",k," bases):\n   either reduce basis.maxdim or collapse categories",sep=""))

  if(basis.maxdim < 1) stop(" basis.maxdim must be greater than or equal to 1")

  console <- newLineConsole()
  console <- printPush("Working...",console = console)

  ## Exhaustive evaluation over all combinations of K, search over
  ## lambda for each combination

  K.mat <- matrix.combn(0:basis.maxdim,num.x)
  nrow.K.mat <- NROW(K.mat)
  cv.min.vec <- numeric(nrow.K.mat)
  basis.vec <- character(nrow.K.mat)
  lambda.mat <- matrix(NA,nrow.K.mat,num.z)

  cv.min <- .Machine$double.xmax

  t2 <- Sys.time() ## placeholder

  output <- list()
  output.restart <- list()

  for(j in 1:nrow.K.mat) {

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
                        degree=degree,
                        nbreak=nbreak,
                        complexity=complexity,
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
                                    nbreak=nbreak,
                                    degree=degree,
                                    complexity=complexity,
                                    basis="additive-tensor")

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      cv.min.vec[j] <- output$value
      basis.vec[j] <- "additive-tensor"

      if(output$value < cv.min) {
        output.opt <- output
        cv.min <- output$value
        K.opt <- K.mat[j,]
        lambda.opt <- output$par
        basis.opt <- "additive-tensor"
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
                        nbreak=nbreak,
                        degree=degree,
                        complexity=complexity,
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
                                    nbreak=nbreak,
                                    degree=degree,
                                    complexity=complexity,
                                    basis="additive")

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      if(output$value < cv.min) {
        output.opt <- output
        cv.min <- output$value
        K.opt <- K.mat[j,]
        lambda.opt <- output$par
        basis.opt <- "additive"
        cv.min.vec[j] <- output$value
        basis.vec[j] <- "additive"
      }

      lambda.mat[j,] <- output$par

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
                        nbreak=nbreak,
                        degree=degree,
                        complexity=complexity,
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
                                    degree=degree,
                                    nbreak=nbreak,
                                    complexity=complexity,
                                    basis="tensor")

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      if(output$value < cv.min) {
        output.opt <- output
        cv.min <- output$value
        K.opt <- K.mat[j,]
        lambda.opt <- output$par
        basis.opt <- "tensor"
        cv.min.vec[j] <- output$value
        basis.vec[j] <- "tensor"
      }

      lambda.mat[j,] <- output$par

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
                        degree=degree,
                        nbreak=nbreak,
                        complexity=complexity,
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
                                    degree=degree,
                                    nbreak=nbreak,
                                    complexity=complexity,
                                    basis=basis)

          }

          if(output.restart$value < output$value) output <- output.restart

        }

      } ## end restarts

      if(output$value < cv.min) {
        output.opt <- output
        cv.min <- output$value
        K.opt <- K.mat[j,]
        lambda.opt <- output$par
        basis.opt <- basis
      }

      basis.vec[j] <- basis
      cv.min.vec[j] <- output$value
      lambda.mat[j,] <- output$par

    }

    t2 <- Sys.time()

  }

  console <- printClear(console)
  console <- printPop(console)

  if(any(K.opt==basis.maxdim)) warning(paste(" optimal K equals search maximum (", basis.maxdim,"): rerun with larger basis.maxdim",sep=""))

  crscv(K=K.opt,
        I=NULL,
        basis=basis.opt,
        basis.vec=basis.vec,
        basis.maxdim=basis.maxdim,
        restarts=restarts,
        K.mat=K.mat,
        lambda=lambda.opt,
        lambda.mat=lambda.mat,
        cv.func=cv.min,
        cv.func.vec=as.matrix(cv.min.vec))

}
