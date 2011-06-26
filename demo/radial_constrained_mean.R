## $Id: radial_constrained_mean.R,v 1.1 2011/06/25 15:09:36 jracine Exp jracine $

## Code to conduct restricted regression splines on evaluation
## data. Presumes continuous regressors, accepts an arbitrary number
## of regressors, and accepts arbitrary derivative restrictions.

rm(list=ls())

## Parameters to be set.

set.seed(42)

n <- 1000
n.eval <- 50

x.min <- -5
x.max <- 5

## These will need to be modified if/when you modify Amat and bvec

lower <- 0
upper <- 0.5

## Load libraries

require(crs)
require(quadprog)

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

x1.seq <- seq(x.min,x.max,length=n.eval)
x2.seq <- seq(x.min,x.max,length=n.eval)
data.eval <- data.frame(y=0,expand.grid(x1=x2.seq,x2=x2.seq))

model.unres <- crs(y~x1+x2,
                   basis="auto",
                   data=data.train,
                   nmulti=5)

summary(model.unres)

## Start from uniform weights equal to 1/n. If constraints are
## non-binding these are optimal.

p <- rep(1/n,n)
Dmat <- diag(1,n,n)
dvec <- as.vector(p)

## If you wish to alter the constraints, you need to modify Amat and
## bvec.

## Generate the estimated model computed for the training data. Note -
## we need to premultiply the weights by n and each column must be
## multiplied by y

B <- model.matrix(model.unres$model.lm)
Aymat.res <- n*t(t(B%*%solve(t(B)%*%B)%*%t(B))*data.train$y)

## Here is Amat

Amat <- t(rbind(rep(1,n),
                Aymat.res,
                -Aymat.res))

rm(Aymat.res)

## Here is bvec

bvec <- c(0,
          (rep(lower,n)-fitted(model.unres)),
          (fitted(model.unres)-rep(upper,n)))

## Solve the quadratic programming problem

QP.output <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)

## No longer needed...

rm(Amat,bvec,Dmat,dvec)

## Get the solution and update the uniform weights

w.hat <- QP.output$solution

p.updated <- p + w.hat

## Now estimate the restricted model

data.trans <- data.frame(y=p.updated*n*data.train$y,data.train[,2:ncol(data.train)])
names(data.trans) <- names(data.train) ## Necessary when there is only 1 regressor
model.res <- crs(y~x1+x2,cv="none",
                 degree=model.unres$degree,
                 segments=model.unres$segments,
                 basis=model.unres$basis,                                  
                 data=data.trans,
                 deriv=1)

## That's it!

## Create a 3D perspective plot of the constrained and unconstrained
## surfaces

fitted.unres <- matrix(predict(model.unres,newdata=data.eval), n.eval, n.eval)
fitted.res <- matrix(predict(model.res,newdata=data.eval), n.eval, n.eval)

ylim <- c(min(fitted(model.unres),max(fitted(model.unres))))

zlim <- c(min(fitted(model.unres)),max(fitted(model.unres)))

par(mfrow=c(1,2))

persp(x1.seq, x2.seq,
      fitted.unres,
      main="Unconstrained Regression Spline",
      col="lightblue",
      ticktype="detailed", 
      ylab="X2",
      xlab="X1",
      zlim=zlim,
      zlab="Conditional Expectation",
      theta=300,
      phi=30)

persp(x1.seq, x2.seq,
      fitted.res,
      main="Constrained Regression Spline",
      sub="0 <= g(x) <= 1/2",
      col="lightblue",
      ticktype="detailed", 
      ylab="X2",
      xlab="X1",
      zlim=zlim,
      zlab="Conditional Expectation",
      theta=300,
      phi=30)


