## This function conducts kernel regression spline
## cross-validation. It takes as input a data.frame xz containing a
## mix of numeric and factor predictors and a vector y. A range of
## arguments can be provided, and one can do search on the bandwidths
## and both the degree and knots ("degree-knots") or the degree
## holding the number of knots (segments+1) constant or the number of
## knots (segments+1) holding the degree constant. Two basis types are
## supported ("additive" or "tensor") and the argument "auto" will
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
                  basis=c("additive","tensor","auto"),
                  cv.func=c("cv.ls","cv.gcv","cv.aic"),
                  degree=degree,
                  segments=segments) {

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)
  cv.func <- match.arg(cv.func)  

  ## First define the cv function to be fed to optim

  kernel.type <- match.arg(kernel.type)

  t1 <- Sys.time()

  cv.objc <- function(input,
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
                      cv.func=cv.func) {

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
                           basis=basis,
                           cv.func=cv.func)

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
  ## tensor in this case, so don't waste time doing both).

  if(num.x==1 & basis == "auto") basis <- "additive"

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

  if(complexity=="degree") {
    K.mat <- cbind(K.mat[,1:num.x],matrix(segments,nrow(K.mat),length(segments),byrow=TRUE))
  } else if(complexity=="knots") {
    K.mat <- cbind(matrix(degree,nrow(K.mat),length(degree),byrow=TRUE),K.mat[,1:num.x]+1) 
  } else if(complexity=="degree-knots") {
    K.mat[,(num.x+1):(2*num.x)] <- K.mat[,(num.x+1):(2*num.x)]+1
  }

  ## Strip off all (except one) rows with degree 0 for all continuous
  ## predictors, only leave first row (avoid redundant computation,
  ## this will be computed only once for all predictors having degree
  ## 0, segments 1).

  degree.zero.rows <- which(rowSums(K.mat[,1:num.x,drop=FALSE])==0)[-1]

  if(length(degree.zero.rows) > 0) K.mat <- K.mat[-degree.zero.rows,,drop=FALSE]

  nrow.K.mat <- NROW(K.mat)

  ## Initialize

  cv.vec <- rep(.Machine$double.xmax,nrow.K.mat)

  basis.vec <- character(nrow.K.mat)
  lambda.mat <- matrix(NA,nrow.K.mat,num.z)

  output <- list()
  output.restart <- list()

  t2 <- Sys.time() ## placeholder

  for(j in 1:nrow.K.mat) {

    if(basis=="auto") {

      ## First, basis=="additive"

      output$convergence <- 42

      while(output$convergence != 0) {

        output <- optim(par=runif(num.z),
                        cv.objc,
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
                        basis="additive",
                        cv.func=cv.func)

      }

      if(restarts > 0) {

        for(r in 1:restarts) {

          output.restart$convergence <- 42

          while(output.restart$convergence != 0) {

            output.restart <- optim(par=runif(num.z),
                                    cv.objc,
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
                                    basis="additive",
                                    cv.func=cv.func)

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
                        cv.objc,
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
                        basis="tensor",
                        cv.func=cv.func)

      }

      if(restarts > 0) {

        for(r in 1:restarts) {

          output.restart$convergence <- 42

          while(output.restart$convergence != 0) {

            output.restart <- optim(par=runif(num.z),
                                    cv.objc,
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
                                    basis="tensor",
                                    cv.func=cv.func)

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

      ## Either basis=="additive" or "tensor"

      output$convergence <- 42

      while(output$convergence != 0) {

        output <- optim(par=runif(num.z),
                        cv.objc,
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
                        basis=basis,
                        cv.func=cv.func)

      }

      if(restarts > 0) {

        for(r in 1:restarts) {

          output.restart$convergence <- 42

          while(output.restart$convergence != 0) {

            output.restart <- optim(par=runif(num.z),
                                    cv.objc,
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
                                    basis=basis,
                                    cv.func=cv.func)

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
  K.opt <- K.mat[ocv.vec,,drop=FALSE][1,]
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
        cv.objc=cv.min,
        cv.objc.vec=as.matrix(cv.vec),
        num.x=num.x,
        cv.func=cv.func)

}
