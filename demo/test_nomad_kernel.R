rm(list=ls())

library(crs)

set.seed(42)

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
                    cv="none",
                    kernel=TRUE)

summary(model.kernel)

## Use initial value starting points

x0<-c(1, 1, 1, 2, 0.5)

model.kernel.x0 <- crs(y~x1+x2+z,
                       degree=c(5,5),
                       lambda=c(0.1),
                       basis.maxdim=6,
                       basis="auto",
                       comlexity="degree-knots",
                       cv="nomad",
                       kernel=TRUE,
                       x0=x0)

summary(model.kernel.x0)

## Multiple initial points - x0 will be the first inital point

model.kernel.multiple <- crs(y~x1+x2+z,
                             degree=c(5,5),
                             lambda=c(0.1),
                             basis.maxdim=6,
                             basis="auto",
                             comlexity="degree-knots",
                             cv="nomad",
                             kernel=TRUE,
                             x0=x0,
                             nmulti=10)

summary(model.kernel.multiple)

## We could compare with exhaustive search, but that takes some time
## as numerical search is conducted for each degree/segment combination.

#model.kernel.multiple <- crs(y~x1+x2+z,
#                             degree=c(5,5),
#                             lambda=c(0.1),
#                             basis.maxdim=6,
#                             basis="auto",
#                             comlexity="degree-knots",
#                             cv="exhaustive",
#                             kernel=TRUE,
#                             x0=x0,
#                             nmulti=10)

