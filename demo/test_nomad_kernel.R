rm(list=ls())

library(crs)

set.seed(123)

## Example - simulated data

n <- 10000
num.eval <- 50
x1 <- runif(n)
x2 <- runif(n)
z <- rbinom(n,1,.5)
dgp <- cos(2*pi*x1)+sin(2*pi*x2)+z
z <- factor(z)
y <- dgp + rnorm(n,sd=.5)

## Estimate a model with specified degree and bandwidth

model.kernel <- crs(y~x1+x2+z,
                    degree=c(5,5),
                    lambda=c(0.1),
                    basis.maxdim=6,
                    basis="additive",
                    comlexity="degree-knots",
                    cv=FALSE,
                    kernel=TRUE,
                    nomad=TRUE)

summary(model.kernel)

## Could you check the error, I guess the problem is that I do not
## output K.mat which you will call it to do the calculation later.

x0<-c(1, 1, 1, 2, 0.5)

model.kernel.x0 <- crs(y~x1+x2+z,
                       degree=c(5,5),
                       lambda=c(0.1),
                       basis.maxdim=6,
                       basis="auto",
                       comlexity="degree-knots",
                       cv=TRUE,
                       kernel=TRUE,
                       nomad=TRUE,
                       x0=x0)

summary(model.kernel.x0)

## Multiple initial points - x0 will be the first inital point

model.kernel.multiple <- crs(y~x1+x2+z,
                             degree=c(5,5),
                             lambda=c(0.1),
                             basis.maxdim=6,
                             basis="auto",
                             comlexity="degree-knots",
                             cv=TRUE,
                             kernel=TRUE,
                             nomad=TRUE,
                             x0=x0,
                             nb_mads_runs=10)

summary(model.kernel.multiple)
