## This demo considers nonparametric instrumental regression in a
## setting with one endogenous regressor, one exogenous regressor, and
## one instrument.

## This illustration was made possible by Samuele Centorrino
## <samuele.centorrino@univ-tlse1.fr>

require(crs)

## Turn off screen I/O for crs()

opts <- list("MAX_BB_EVAL"=10000,
             "EPSILON"=.Machine$double.eps,
             "INITIAL_MESH_SIZE"="r1.0e-01",
             "MIN_MESH_SIZE"=sqrt(.Machine$double.eps),
             "MIN_FRAME_SIZE"=sqrt(.Machine$double.eps),
             "DISPLAY_DEGREE"=0)

set.seed(42)

.crs_demo_numeric <- function(prompt, default, env = NULL, choices = NULL) {
  value <- if(!is.null(env)) Sys.getenv(env, unset = "") else ""
  if(!nzchar(value) && interactive()) value <- readline(prompt = prompt)
  if(!nzchar(value)) value <- as.character(default)
  value <- suppressWarnings(as.numeric(value))
  if(!is.finite(value) || (!is.null(choices) && !(value %in% choices))) {
    value <- default
  }
  value
}

## Interactively request number of observations, the method, whether
## to do NOMAD or exhaustive search, and if NOMAD the number of
## multistarts

n <- .crs_demo_numeric("Input the number of observations desired: ",
                       100,
                       "CRS_DEMO_N")
method <- .crs_demo_numeric("Input the method (0=Landweber-Fridman, 1=Tikhonov): ",
                            0,
                            "CRS_DEMO_METHOD",
                            choices = c(0, 1))
method <- ifelse(method==0,"Landweber-Fridman","Tikhonov")
cv <- .crs_demo_numeric("Input the cv method (0=nomad, 1=exhaustive): ",
                        1,
                        "CRS_DEMO_CV",
                        choices = c(0, 1))
cv <- ifelse(cv==0,"nomad","exhaustive")
nmulti <- 1
if(cv=="nomad") {
  nmulti <- .crs_demo_numeric("Input the number of multistarts desired (e.g. 10): ",
                              1,
                              "CRS_DEMO_NMULTI")
}
num.eval <- .crs_demo_numeric("Input the number of evaluation observations desired (e.g. 50): ",
                              25,
                              "CRS_DEMO_NUM_EVAL")

v  <- rnorm(n,mean=0,sd=.27)
eps <- rnorm(n,mean=0,sd=0.05)
u <- -0.5*v + eps
w <- rnorm(n,mean=0,sd=1)
z <- 0.2*w + v
x <- rnorm(n)

## In Darolles et al (2011) there exist two DGPs. The first is
## phi(z)=z^2. Here we add an exogenous regressor.

phi <- function(z) { z^2 }
eyz <- function(z) { z^2 -0.325*z }

y <- phi(z) + 0.05*x^3 + u

## In evaluation data sort z for plotting and hold x constant at its
## median

## Setting cv.threshold = 0 forces NOMAD search instead of exhaustive search
## when no categorical predictors are present. This avoids unnecessary
## evaluation of all degree/segment combinations in the examples and, for
## crsiv() and crsivderiv(), ensures that the warm-start strategy is used.
model.iv <- crsiv(y=y,z=z,w=w,x=x,cv=cv,nmulti=nmulti,method=method,deriv=1,cv.threshold=0)

## Setting cv.threshold = 0 forces NOMAD search instead of exhaustive search
## when no categorical predictors are present. This avoids unnecessary
## evaluation of all degree/segment combinations in the examples and, for
## crsiv() and crsivderiv(), ensures that the warm-start strategy is used.
model.noniv <- crs(y~z+x,cv=cv,nmulti=nmulti,deriv=1,opts=opts,cv.threshold=0)

summary(model.iv)

# Perspective plot
z.seq <- seq(min(z),max(z),length=num.eval)
x.seq <- seq(min(x),max(x),length=num.eval)
x.grid <- expand.grid(z.seq,x.seq)
newdata <- data.frame(z=x.grid[,1],x=x.grid[,2])
z.iv <- matrix(predict(model.iv,newdata=newdata),num.eval,num.eval)
z.noniv <- matrix(predict(model.noniv,newdata=newdata),num.eval,num.eval)
zlim <- c(min(z.iv,z.noniv),max(z.iv,z.noniv))
persp(x=z.seq,y=x.seq,z=z.iv,
      xlab="Z",ylab="X",zlab="Y",
      zlim=zlim,
      ticktype="detailed",
      col=FALSE,
      border="red",
      main="phi(z,x)",
      theta=45,phi=45)
par(new=TRUE)
persp(x=z.seq,y=x.seq,z=z.noniv,
      xlab="Z",ylab="X",zlab="Y",
      zlim=zlim,
      ticktype="detailed",
      col=FALSE,
      border="blue",
      main="phi(z,x)",
      theta=45,phi=45)

## Perspective plot - derivative wrt z
z.iv <- matrix(attr(predict(model.iv,newdata=newdata),"deriv.mat")[,1],num.eval,num.eval)
z.noniv <- matrix(attr(predict(model.noniv,newdata=newdata),"deriv.mat")[,1],num.eval,num.eval)
zlim <- c(min(z.iv,z.noniv),max(z.iv,z.noniv))
persp(x=z.seq,y=x.seq,z=z.iv,
      xlab="Z",ylab="X",zlab="Y",
      zlim=zlim,
      ticktype="detailed",      
      border="red",
      col=FALSE,
      main="d g(z,x)/d z (x=med(x))",
      theta=45,phi=45)
par(new=TRUE)
persp(x=z.seq,y=x.seq,z=z.noniv,
      xlab="Z",ylab="X",zlab="Y",
      zlim=zlim,
      ticktype="detailed",      
      border="blue",
      col=FALSE,
      main="d g(z,x)/d z (x=med(x))",
      theta=45,phi=45)

## Perspective plot - derivative wrt x
z.iv <- matrix(attr(predict(model.iv,newdata=newdata),"deriv.mat")[,2],num.eval,num.eval)
z.noniv <- matrix(attr(predict(model.noniv,newdata=newdata),"deriv.mat")[,2],num.eval,num.eval)
zlim <- c(min(z.iv,z.noniv),max(z.iv,z.noniv))
persp(x=z.seq,y=x.seq,z=z.iv,
      xlab="Z",ylab="X",zlab="Y",
      zlim=zlim,
      ticktype="detailed",      
      border="red",
      col=FALSE,
      main="d g(z,x)/d x (z=med(z))",
      theta=45,phi=45)
par(new=TRUE)
persp(x=z.seq,y=x.seq,z=z.noniv,
      xlab="Z",ylab="X",zlab="Y",
      zlim=zlim,
      ticktype="detailed",      
      border="blue",
      col=FALSE,
      main="d g(z,x)/d x (z=med(z))",
      theta=45,phi=45)
