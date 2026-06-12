## This illustration considers the `radial function' and plots the
## results in a 3D perspective plot.

require(crs)

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

## Interactively request number of observations, whether to do NOMAD
## or exhaustive search, and if NOMAD the number of multistarts

n <- .crs_demo_numeric("Input the number of observations desired: ",
                       100,
                       "CRS_DEMO_N")
cv <- .crs_demo_numeric("Input the cv method (0=nomad, 1=exhaustive): ",
                        0,
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

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

model <- crs(y~x1+x2,
             cv=cv,
             complexity="degree-knots",
             knots="uniform",
             deriv=1,
             cv.func="cv.aic",
             nmulti=nmulti)

summary(model)

## Perspective plot
plot(model,perspective=TRUE,view="fixed",main="Conditional Mean")

## Derivative surface plots
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])

## Perspective plot - derivative wrt x1
z <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,1],num.eval,num.eval)

persp(x=x1.seq,y=x2.seq,z=z,
      xlab="X1",ylab="X2",zlab="Y",
      ticktype="detailed",
      col=grDevices::adjustcolor("red",alpha.f=0.35),
      border=grDevices::adjustcolor("red",alpha.f=0.60),
      main="d g(x1,x2)/d x1 (x2=med(x2))",
      theta=45,phi=45)

## Perspective plot - derivative wrt x2
z <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,2],num.eval,num.eval)

persp(x=x1.seq,y=x2.seq,z=z,
      xlab="X1",ylab="X2",zlab="Y",
      ticktype="detailed",
      col=grDevices::adjustcolor("red",alpha.f=0.35),
      border=grDevices::adjustcolor("red",alpha.f=0.60),
      main="d g(x1,x2)/d x2 (x1=med(x1))",
      theta=45,phi=45)
