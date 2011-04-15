## This presents some 3D plots

rm(list=ls())

library(crs)

set.seed(42)

n <- 10000

num.eval <- 50
x1 <- runif(n)
x2 <- runif(n)
z <- rbinom(n,1,.5)
dgp <- cos(2*pi*x1)+sin(2*pi*x2)+z
z <- factor(z)
y <- dgp + rnorm(n,sd=.5)

## Estimate a model with pre-specified degree and bandwidth

model <- crs(y~x1+x2+z,
             degree=c(3,4),
             segments=c(3,2),
             lambda=c(0.1),
             basis="additive",
             cv=FALSE,
             kernel=TRUE,
             knots="uniform",
             deriv=1)

summary(model)

## Perspective plot
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(0,num.eval**2),levels=c(0,1)))
z0 <- matrix(predict(model,newdata=newdata),num.eval,num.eval)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(1,num.eval),levels=c(0,1)))
z1 <- matrix(predict(model,newdata=newdata),num.eval,num.eval)
zlim=c(min(z0,z1),max(z0,z1))
persp(x=x1.seq,y=x2.seq,z=z0,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      ticktype="detailed",      
      border="red",
      main="Conditional Mean",
      theta=45,phi=45)
par(new=TRUE)
persp(x=x1.seq,y=x2.seq,z=z1,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      theta=45,phi=45,
      ticktype="detailed",
      border="blue")

dev.new()

## Perspective plot - derivative wrt x1
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(0,num.eval**2),levels=c(0,1)))
z0 <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,1],num.eval,num.eval)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(1,num.eval),levels=c(0,1)))
z1 <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,1],num.eval,num.eval)
zlim=c(min(z0,z1),max(z0,z1))
persp(x=x1.seq,y=x2.seq,z=z0,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      ticktype="detailed",      
      border="red",
      main="d g(x1,x2)/dx1 (x2=med(x2))",
      theta=45,phi=45)
par(new=TRUE)
persp(x=x1.seq,y=x2.seq,z=z1,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      theta=45,phi=45,
      ticktype="detailed",
      border="blue")

dev.new()

## Perspective plot - derivative wrt x2
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(0,num.eval**2),levels=c(0,1)))
z0 <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,2],num.eval,num.eval)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(1,num.eval),levels=c(0,1)))
z1 <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,2],num.eval,num.eval)
zlim=c(min(z0,z1),max(z0,z1))
persp(x=x1.seq,y=x2.seq,z=z0,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      ticktype="detailed",      
      border="red",
      main="d g(x1,x2)/dx2 (x1=med(x1))",
      theta=45,phi=45)
par(new=TRUE)
persp(x=x1.seq,y=x2.seq,z=z1,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      theta=45,phi=45,
      ticktype="detailed",
      border="blue")

