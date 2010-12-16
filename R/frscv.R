## This function conducts factor regression spline
## cross-validation. It takes as input a data.frame xz containing a
## mix of numeric and factor predictors and a vector y. A range of
## arguments can be provided, and one can do search on both the degree
## and knots ("degree-knots") or the degree holding the number of
## knots (segments+1) constant or the number of knots (segments+1)
## holding the degree constant. A variety of basis types are supported
## (functional anova "additive-tensor", "additive", or "tensor") and
## the argument "auto" will evaluate choose the basis type
## automatically.

## Currently search is exhaustive taking basis.maxdim as the maximum
## number of the spline degree (0,1,...) and number of segments
## (1,2,...). This is a quadratic integer programming problem so
## ideally I require an IQP (MIQP for kernel-weighting)
## solver. Currently in R there is no such beast.

frscv <- function(xz,
                  y,
                  basis.maxdim=5,
                  complexity=c("degree-knots","degree","knots"),
                  knots=c("quantiles","uniform"),
                  basis=c("additive-tensor","additive","tensor","auto"),
                  cv.norm=c("L2","L1"),
                  degree=degree,
                  segments=segments) {

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)
  cv.norm <- match.arg(cv.norm)

  t1 <- Sys.time()

  cv.func <- function(input,
                      x,
                      z=NULL,
                      y,
                      restart,
                      num.restarts,
                      basis.maxdim,
                      j=NULL,
                      nrow.KI.mat=NULL,
                      t2=NULL,
                      complexity=complexity,
                      knots=knots,
                      basis=basis,
                      cv.norm=cv.norm) {

    if(missing(input) || missing(x) || missing(y) || missing(basis.maxdim)) stop(" you must provide input, x, y, and basis.maxdim")

    ## Presumes x (continuous predictors) exist, but z
    ## (ordinal/nominal factors) can be optional

    n <- length(y)
    num.x <- NCOL(x)
    
    if(NROW(y) != NROW(x)) stop(" x and y have differing numbers of observations")

    ## K is a matrix, column 1 degree, column 2 segments, either or
    ## both can be determined via cv so need to take care to allow
    ## user to select knots (degree fixed), degree (knots fixed), or
    ## both degree and knots. The values used to evaluate the cv
    ## function are passed below.
    
    K <- round(cbind(input[1:num.x],input[(num.x+1):(2*num.x)]))

    if(!is.null(z)) {
      num.z <- NCOL(z)
      I <- round(input[(2*num.x+1):(2*num.x+num.z)])
    } else {
      num.z <- 0
      I <- NULL
    }

    cv <- cv.factor.spline(x=x,
                           y=y,
                           z=z,
                           K=K,
                           I=I,
                           knots=knots,
                           basis=basis)

    ## Some i/o unless options(crs.messages=FALSE)

    ## Degree is first column of K K[,1], segments second column K[,2]
    ## - could create a tmp vector for i/o, or could switch

    console <<- printClear(console)

    ## Format function...
    fw.format.2 <- function(input) sapply(input,sprintf,fmt="%#.2f")

    ## i/o depends on whether we are cross-validating degree, knots,
    ## or both

    if(complexity=="degree") {
      if(!is.null(j)) {
        if(j==1) {
          tmp.1 <- paste(j,"/",nrow.KI.mat,", d[1]=",K[1,1],sep="")
        } else {
          dt <- (t2-t1)*(nrow.KI.mat-j+1)/j
          tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                         fw.format.2(as.numeric((t2-t1),units="mins")),
                         "m",sep="")
          tmp.1 <- paste(j,"/",nrow.KI.mat,tmp.0,", d[1]=",K[1,1],sep="")
        }
      } else {
        tmp.1 <- paste("d[1]=", K[1,1],sep="")
      }
      if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", d[", i, "]=", K[i,1],sep="")
    } else if(complexity=="knots") {
      if(!is.null(j)) {
        if(j==1) {
          tmp.1 <- paste(j,"/",nrow.KI.mat,", s[1]=",K[1,2],sep="")
        } else {
          dt <- (t2-t1)*(nrow.KI.mat-j+1)/j
          tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                         fw.format.2(as.numeric((t2-t1),units="mins")),
                         "m",sep="")
          tmp.1 <- paste(j,"/",nrow.KI.mat,tmp.0,", s[1]=",K[1,2],sep="")
        }
      } else {
        tmp.1 <- paste("s[1]=", K[1,2],sep="")
      }
      if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", s[", i, "]=", K[i,2],sep="")
    } else if(complexity=="degree-knots") {
      if(!is.null(j)) {
        if(j==1) {
          tmp.1 <- paste(j,"/",nrow.KI.mat,", d[1]=",K[1,1],sep="")
        } else {
          dt <- (t2-t1)*(nrow.KI.mat-j+1)/j
          tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                         fw.format.2(as.numeric((t2-t1),units="mins")),
                         "m",sep="")
          tmp.1 <- paste(j,"/",nrow.KI.mat,tmp.0,", d[1]=",K[1,1],sep="")
        }
      } else {
        tmp.1 <- paste("k[1]=", K[1,1],sep="")
      }
      if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", d[", i, "]=", K[i,1],sep="")
      for(i in 1:num.x) tmp.1 <- paste(tmp.1, ", s[", i, "]=", K[i,2],sep="")      
    }

    ## For i/o for z variables...
    
    if(num.z > 0) for(i in 1:num.z) tmp.1 <- paste(tmp.1, ", I[", i, "]=", I[i],sep="")
    tmp.3 <- paste(", cv=", format(cv,digits=6), sep="")
    if(num.restarts > 0) {
      tmp.2 <- paste(", rs=", restart, "/", num.restarts,sep="")
      msg <- paste(tmp.1,tmp.2,tmp.3,sep="")
    } else {
      msg <- paste(tmp.1,tmp.3,sep="")
    }

    console <<- printPush(msg,console = console)

    return(cv)

  }

  console <- newLineConsole()
  console <- printPush("Working...",console = console)

  ## Take data frame x and parse into factors (z) and numeric (x)

  if(!is.data.frame(xz)) stop(" xz must be a data frame")

  xztmp <- splitFrame(xz)
  x <- xztmp$x
  z <- xztmp$z
  if(is.null(z)) {
    include <- NULL
    num.z <- 0
  } else {
    num.z <- NCOL(z)
  }

  num.x <- ncol(x)
  n <- nrow(x)

  if(missing(x) || missing(y)) stop (" you must provide x and y")

  if(complexity=="degree") {
    if(missing(segments)) stop("segments missing for cross-validation of spline degree")
    if(length(segments)!=num.x) stop(" segments vector must be the same length as x")  
  } else if(complexity=="knots") {
    if(missing(degree)) stop("degree missing for cross-validation of number of spline knots")
    if(length(degree)!=num.x) stop(" degree vector must be the same length as x")
  }

  ## For factor regression spline, if there is only one predictor
  ## (i.e. num.x + num.z = 1) disable auto, set to additive (which is
  ## additive-tensor and tensor in this case, so don't waste time
  ## doing all three).

  if(num.x+num.z==1 & basis == "auto") basis <- "additive"

  ## This call will trap zero/negative degrees of freedom immediately
  ## rather than letting cv proceed only to be halted after the lower
  ## dimension models have been estimated when this occurs.

  if(complexity=="degree") {
    if(basis == "auto") {
      k <- max(c(ncol(prod.spline(x=x,z=z,K=cbind(rep(basis.maxdim,num.x),segments),I=rep(1,num.z),knots=knots,basis="additive-tensor")),
                 ncol(prod.spline(x=x,z=z,K=cbind(rep(basis.maxdim,num.x),segments),I=rep(1,num.z),knots=knots,basis="additive")),
                 ncol(prod.spline(x=x,z=z,K=cbind(rep(basis.maxdim,num.x),segments),I=rep(1,num.z),knots=knots,basis="tensor"))))
    } else {
      k <- ncol(prod.spline(x=x,z=z,K=cbind(rep(basis.maxdim,num.x),segments),I=rep(1,num.z),knots=knots,basis=basis))
    }
  } else if(complexity=="knots") {
    if(basis == "auto") {
      k <- max(c(ncol(prod.spline(x=x,z=z,K=cbind(degree,rep(basis.maxdim,num.x)),I=rep(1,num.z),knots=knots,basis="additive-tensor")),
                 ncol(prod.spline(x=x,z=z,K=cbind(degree,rep(basis.maxdim,num.x)),I=rep(1,num.z),knots=knots,basis="additive")),
                 ncol(prod.spline(x=x,z=z,K=cbind(degree,rep(basis.maxdim,num.x)),I=rep(1,num.z),knots=knots,basis="tensor"))))
    } else {
      k <- ncol(prod.spline(x=x,z=z,K=cbind(degree,rep(basis.maxdim,num.x)),I=rep(1,num.z),knots=knots,basis=basis))
    }
  } else if(complexity=="degree-knots"){
    if(basis == "auto") {
      k <- max(c(ncol(prod.spline(x=x,z=z,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),I=rep(1,num.z),knots=knots,basis="additive-tensor")),
                 ncol(prod.spline(x=x,z=z,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),I=rep(1,num.z),knots=knots,basis="additive")),
                 ncol(prod.spline(x=x,z=z,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),I=rep(1,num.z),knots=knots,basis="tensor"))))
    } else {
      k <- ncol(prod.spline(x=x,z=z,K=matrix(2*rep(basis.maxdim,num.x),num.x,2),I=rep(1,num.z),knots=knots,basis=basis))
    }
  }

  df <- n - k

  if(df <= 0) {
    stop(paste(" maximum basis dimension (",k,") would equal/exceed sample size (",n,")\n   perhaps use basis=\"additive\" or else decrease basis.maxdim",sep=""))
  } else if(df <= 10) {
    warning(paste(" maximum basis dimension (",k,") and sample size (",n,") close",sep=""))
  }

  if(basis.maxdim < 1) stop(" basis.maxdim must be greater than or equal to 1")

  ## Need to append 0,1 for I (in, out). XXX Pad row to left with
  ## degree vector or in between x and z with segment vector then
  ## input remains a vector XXX need to do... this ought to be most
  ## transparent and the remaining code will not need to be changed XXX

  if(complexity!="degree-knots") {
    if(!is.null(z)) {
      KI.mat <- matrix.combn(0:basis.maxdim,num.x,num.z)
    } else {
      KI.mat <- matrix.combn(0:basis.maxdim,num.x)
    }
  } else {
    if(!is.null(z)) {
      KI.mat <- matrix.combn(0:basis.maxdim,2*num.x,num.z)
    } else {
      KI.mat <- matrix.combn(0:basis.maxdim,2*num.x)
    }
  }
  
  nrow.KI.mat <- NROW(KI.mat)
  basis.vec <- character(nrow.KI.mat)
  cv.vec <- numeric(nrow.KI.mat)

  for(j in 1:nrow.KI.mat) {

    ## Initialize    

    cv.vec[j] <- .Machine$double.xmax

    if(complexity=="degree") {
      if(num.z==0) {
        input.j <- c(KI.mat[j,1:num.x],segments)
      } else {
        input.j <- c(KI.mat[j,1:num.x],segments,KI.mat[j,(num.x+1):(num.x+num.z)])
      }
    } else if(complexity=="knots") {
      if(num.z==0) {
        input.j <- c(degree,KI.mat[j,1:num.x]+1)
      } else {
        input.j <- c(degree,KI.mat[j,1:num.x]+1,KI.mat[j,(num.x+1):(num.x+num.z)])
      }
    } else if(complexity=="degree-knots") {
      KI.mat[j,(num.x+1):(2*num.x)] <- KI.mat[j,(num.x+1):(2*num.x)]+1
      input.j <- KI.mat[j,]
    }

    if(basis=="auto") {

      output <- cv.func(input=input.j,
                        x=x,
                        y=y,
                        z=z,
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis="additive-tensor")

      if(output < cv.vec[j]) {
        cv.vec[j] <- output
        basis.vec[j] <- "additive-tensor"
      }


      output <- cv.func(input=input.j,
                        x=x,
                        y=y,
                        z=z,
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis="additive")

      if(output < cv.vec[j]) {
        cv.vec[j] <- output
        basis.vec[j] <- "additive"
      }        

      output <- cv.func(input=input.j,
                        x=x,
                        y=y,
                        z=z,
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis="tensor")

      if(output < cv.vec[j]) {
        cv.vec[j] <- output
        basis.vec[j] <- "tensor"
      }

    } else {

      ## not auto, so use either "additive-tensor" or "additive" or "tensor"

      output <- cv.func(input=input.j,
                        x=x,
                        y=y,
                        z=z,
                        basis.maxdim=basis.maxdim,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis=basis)

      if(output < cv.vec[j]) {
        cv.vec[j] <- output
        basis.vec[j] <- basis
      }

    }

  }

  ## Sort on cv.vec

  ocv.vec <- order(cv.vec)

  cv.min <- cv.vec[ocv.vec][1]
  if(complexity=="degree") {
    if(num.z==0) {
      K.opt <- c(KI.mat[ocv.vec,1:num.x,drop=FALSE][1,],segments)
    } else {
      K.opt <- c(KI.mat[ocv.vec,1:num.x,drop=FALSE][1,],segments,KI.mat[ocv.vec,(num.x+1):(num.x+num.z),drop=FALSE][1,])
    }
  } else if(complexity=="knots") {
    if(num.z==0) {
      K.opt <- c(degree,KI.mat[ocv.vec,1:num.x,drop=FALSE][1,]+1)
    } else {
      K.opt <- c(degree,KI.mat[ocv.vec,1:num.x,drop=FALSE][1,]+1,KI.mat[ocv.vec,(num.x+1):(num.x+num.z),drop=FALSE][1,])
    }
  } else if(complexity=="degree-knots") {
    K.opt <- KI.mat[ocv.vec,,drop=FALSE][1,]
  }
  basis.opt <- basis.vec[ocv.vec][1]
  degree <- K.opt[1:num.x]
  segments <- K.opt[(num.x+1):(2*num.x)]

  if(!is.null(z)) I.opt <- K.opt[(2*num.x+1):(2*num.x+num.z)]

  console <- printClear(console)
  console <- printPop(console)

  if(any(K.opt==basis.maxdim)) warning(paste(" optimal K equals search maximum (", basis.maxdim,"): rerun with larger basis.maxdim",sep=""))

  if(is.null(z)) I.opt <- NULL

  crscv(K=K.opt,
        I=I.opt,
        basis=basis.opt,
        basis.vec=basis.vec,
        basis.maxdim=basis.maxdim,
        complexity=complexity,
        knots=knots,
        degree=degree,
        segments=segments,
        K.mat=KI.mat,
        restarts=NULL,
        lambda=NULL,
        lambda.mat=NULL,
        cv.func=cv.min,
        cv.func.vec=as.matrix(cv.vec))

}
