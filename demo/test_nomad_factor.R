rm(list=ls())

library(crs)

set.seed(42)

n <- 1000
num.eval <- 50

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

## no initial points, generate by sample

model_nomad <- crs(y~x1+x2,
                   basis="auto",
                   cv=TRUE,
                   complexity="degree-knots",
                   basis.maxdim=5,
                   knots="uniform",
                   deriv=1,
                   cv.func="cv.aic",
                   nomad=TRUE)

summary(model_nomad)

## Initial points, please note that the length of x0 should be the
## same as the number of variables in nomad. This will depend on the
## parameter complexity and the number of predictors.

x0<-c(3, 3, 2, 4)

## Single inital point

model_nomad <- crs(y~x1+x2,
                   basis="auto",
                   cv=TRUE,
                   complexity="degree-knots",
                   basis.maxdim=5,
                   knots="uniform",
                   deriv=1,
                   cv.func="cv.aic", 
                   nomad=TRUE, 
                   x0=x0)

summary(model_nomad)

## Multiple initial points, x0 is not provided, if best_x.txt exists,
## it will be read in as the first inital point, otherwise, points
## will be generated randomly.

model_multi_nomad <- crs(y~x1+x2,
                         basis="auto",
                         cv=TRUE,
                         complexity="degree-knots",
                         basis.maxdim=5,
                         knots="uniform",
                         deriv=1,
                         cv.func="cv.aic", 
                         nomad=TRUE,
                         x0 = x0, 
                         nb_mads_runs=10)

summary(model_multi_nomad)

## Single inital point

model <- crs(y~x1+x2,
             basis="auto",
             cv=TRUE,
             complexity="degree-knots",
             basis.maxdim=5,
             knots="uniform",
             deriv=1,
             cv.func="cv.aic",
             nomad=FALSE)

summary(model)
