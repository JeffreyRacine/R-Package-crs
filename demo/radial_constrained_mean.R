## Code to conduct restricted regression splines on evaluation
## data. Presumes continuous regressors, accepts an arbitrary number
## of regressors, and accepts arbitrary derivative restrictions.

## IMPORTANT NOTE - the code that follows is only valid for the tensor
## basis (basis="tensor")

## Load libraries

require(crs)
require(quadprog)

## Parameters to be set.

set.seed(42)

n <- 1000

x.min <- -5
x.max <- 5

## These will need to be modified if/when you modify Amat and bvec

lower <- 0
upper <- 0.5

## IMPORTANT - you must be careful to NOT read data from environment -
## this appears to work - create a data frame.

## IMPORTANT - code that follows presumes y is the first variable in
## the data frame and all remaining variables are regressors used for
## the estimation.

## Generate a DGP, or read in your own data and create y, x1,
## etc. When you change this by adding or removing variables you need
## to change `data', `rm(...)', and `bw <- ...'. After that all code
## will need no modification.

x1 <- runif(n,x.min,x.max)
x2 <- runif(n,x.min,x.max)

y <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2) + rnorm(n,sd=.1)

data.train <- data.frame(y,x1,x2)

rm(y,x1,x2)

model.unres <- crs(y~x1+x2,
                   data=data.train,
                   nmulti=5,
                   basis="tensor")

summary(model.unres)

## If you wish to alter the constraints, you need to modify Amat and
## bvec.

## Construct the response-scaled constraint operator.

Aymat.res <- crshat(model.unres, y=data.train$y, output="constraint")

## Here is Amat

Amat <- cbind(Aymat.res,
              -Aymat.res)

rm(Aymat.res)

## Here is bvec

bvec <- c(rep(lower,n),
          -rep(upper,n))

## Solve the quadratic programming problem

QP.output <- solve.QP(Dmat=diag(n),dvec=rep(1,n),Amat=Amat,bvec=bvec)

if(is.nan(QP.output$value)) stop(" solve.QP failed. Try smoother curve (larger bandwidths or polynomial order)")

## No longer needed...

rm(Amat,bvec)

## Get the solution

p.hat <- QP.output$solution

## Now estimate the restricted model

data.trans <- data.frame(y=p.hat*data.train$y,data.train[,2:ncol(data.train),drop=FALSE])

model.res <- crs(y~x1+x2,cv="none",
                 degree=model.unres$degree,
                 segments=model.unres$segments,
                 basis=model.unres$basis,
                 data=data.trans,
                 deriv=1)

## That's it!

## Create perspective plots of the constrained and unconstrained surfaces.

par(mfrow=c(1,2))

plot(model.unres,perspective=TRUE,view="fixed",
     main="Unconstrained Regression Spline",
     zlab="Conditional Expectation")

plot(model.res,perspective=TRUE,view="fixed",
     main="Constrained Regression Spline",
     sub="0 <= g(x1,x2) <= 1/2",
     zlab="Conditional Expectation")

par(mfrow=c(1,1))
