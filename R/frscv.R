frscv <- function(xz,
                  y,
                  max.K=10,
                  basis = c("additive-tensor","additive","tensor","auto"),
                  cv.norm=c("L2","L1")) {

  basis <- match.arg(basis)
  cv.norm <- match.arg(cv.norm)  

  t1 <- Sys.time()

  cv.func <- function(input,
                      x,
                      z=NULL,
                      y,
                      restart,
                      num.restarts,
                      max.K,
                      j=NULL,
                      nrow.KI.mat=NULL,
                      t2=NULL,
                      basis=basis,
                      cv.norm=cv.norm) {

    if(missing(input) || missing(x) || missing(y) || missing(max.K)) stop(" you must provide input, x, y, and max.K")

    ## Presumes x (continuous predictors) exist, but z
    ## (ordinal/nominal factors) can be optional

    n <- length(y)
    num.x <- NCOL(x)
    if(NROW(y) != NROW(x)) stop(" x and y have differing numbers of observations")
    K <- round(input[1:num.x])
    if(!is.null(z)) {
      num.z <- NCOL(z)
      I <- round(input[(num.x+1):(num.x+num.z)])
    } else {
      num.z <- 0
      I <- NULL
    }

    cv <- cv.factor.spline(x=x,
                           y=y,
                           z=z,
                           K=K,
                           I=I,
                           basis=basis)

    ## Some i/o unless options(crs.messages=FALSE)

    console <<- printClear(console)

    fw.format.2 <- function(input) sapply(input,sprintf,fmt="%#.2f")
    if(!is.null(j)) {
      if(j==1) {
        tmp.1 <- paste(j,"/",nrow.KI.mat,", k[1]=",K[1],sep="")
      } else {
        dt <- (t2-t1)*(nrow.KI.mat-j+1)/j
        tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                       fw.format.2(as.numeric((t2-t1),units="mins")),
                       "m",sep="")
        tmp.1 <- paste(j,"/",nrow.KI.mat,tmp.0,", k[1]=",K[1],sep="")
      }
    } else {
      tmp.1 <- paste("k[1]=", K[1],sep="")
    }
    if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", k[", i, "]=", K[i],sep="")
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
    num.z <- NULL
  } else {
    num.z <- NCOL(z)
  }

  num.x <- ncol(x)
  n <- nrow(x)

  if(missing(x) || missing(y)) stop (" you must provide x and y")

  ## This call will trap an error immediately rather than letting cv
  ## proceed only to be halted when this occurs

  k <-  max.K*num.x
  if(!is.null(num.z)) {
    num.bases.z <- numeric(length=num.z)
    for(i in 1:num.z) num.bases.z[i] <- (length(unique(as.numeric(z[,i])))-1)
    k <- k + sum(num.bases.z)    
    k <- k + ifelse(basis!="additive",max.K^num.x*prod(num.bases.z),0)
  } else {
    k <- k + ifelse(basis!="additive",max.K^num.x,0)
  }

  df <- n - k

  if(df <= 0) {
    stop(paste(" maximum spline dimension (",k,") equals/exceeds sample size (",n,")\n   perhaps use basis=\"additive\" or else decrease degree.max",sep=""))
  } else if(df <= 10) {
    warning(paste(" maximum spline dimension (",k,") and sample size (",n,") close",sep=""))
  }

  if(max.K < 1) stop(" max.K must be greater than or equal to 1")

  ## Need to append 0,1 for I (in, out)
  if(!is.null(z)) {
    KI.mat <- matrix.combn(0:max.K,num.x,num.z)
  } else {
    KI.mat <- matrix.combn(0:max.K,num.x)
  }
  nrow.KI.mat <- NROW(KI.mat)
  basis.vec <- character(nrow.KI.mat)
  cv.min.vec <- numeric(nrow.KI.mat)

  cv.min <- .Machine$double.xmax

  for(j in 1:nrow.KI.mat) {

    if(basis=="auto") {

      output <- cv.func(input=KI.mat[j,],
                        x=x,
                        y=y,
                        z=z,
                        max.K=max.K,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        basis="additive-tensor")

      if(output < cv.min) {
        cv.min <- output
        K.opt <- KI.mat[j,1:num.x]
        basis.opt <- "additive-tensor"
        basis.vec[j] <- "additive-tensor"
        if(!is.null(z)) I.opt <- KI.mat[j,(num.x+1):(num.x+num.z)]
      }

      output <- cv.func(input=KI.mat[j,],
                        x=x,
                        y=y,
                        z=z,
                        max.K=max.K,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        basis="additive")

      if(output < cv.min) {
        cv.min <- output
        K.opt <- KI.mat[j,1:num.x]
        basis.opt <- "additive"
        basis.vec[j] <- "additive"
        if(!is.null(z)) I.opt <- KI.mat[j,(num.x+1):(num.x+num.z)]
      }

    } else {

      ## not auto, so use either "additive-tensor" or "additive"

      output <- cv.func(input=KI.mat[j,],
                        x=x,
                        y=y,
                        z=z,
                        max.K=max.K,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        basis=basis)

      if(output < cv.min) {
        cv.min <- output
        K.opt <- KI.mat[j,1:num.x]
        basis.opt <- basis
        basis.vec[j] <- basis
        if(!is.null(z)) I.opt <- KI.mat[j,(num.x+1):(num.x+num.z)]
      }

    }

    cv.min.vec[j] <- output

  }

  console <- printClear(console)
  console <- printPop(console)

  if(any(K.opt==max.K)) warning(paste(" optimal K equals search maximum (", max.K,"): rerun with larger max.K",sep=""))

  if(is.null(z)) I.opt <- NULL

  crscv(K=K.opt,
        I=I.opt,
        basis=basis.opt,
        basis.vec=basis.vec,
        max.K=max.K,
        K.mat=KI.mat,
        restarts=NULL,
        lambda=NULL,
        lambda.mat=NULL,
        cv.func=cv.min,
        cv.func.vec=as.matrix(cv.min.vec))

}
