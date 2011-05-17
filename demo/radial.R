rm(list=ls())

library(crs)

set.seed(42)

n <- 1000
num.eval <- 50

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

model <- crs(y~x1+x2,
             basis="auto",
             cv=TRUE,
             complexity="degree-knots",
             basis.maxdim=5,
             knots="uniform",
             deriv=1,
             cv.func="cv.aic")

summary(model)

# Perspective plot
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])
z <- matrix(predict(model,newdata=newdata),num.eval,num.eval)
persp(x=x1.seq,y=x2.seq,z=z,
      xlab="x1",ylab="x2",zlab="y",
      ticktype="detailed",      
      border="red",
      main="Conditional Mean",
      theta=45,phi=45)

## Perspective plot - derivative wrt x1
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])
z <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,1],num.eval,num.eval)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])

persp(x=x1.seq,y=x2.seq,z=z,
      xlab="x1",ylab="x2",zlab="y",
      ticktype="detailed",      
      border="red",
      main="d g(x1,x2)/d x1 (x2=med(x2))",
      theta=45,phi=45)

## Perspective plot - derivative wrt x2
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])
z <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,2],num.eval,num.eval)

persp(x=x1.seq,y=x2.seq,z=z,
      xlab="x1",ylab="x2",zlab="y",
      ticktype="detailed",      
      border="red",
      main="d g(x1,x2)/d x2 (x1=med(x1))",
      theta=45,phi=45)
