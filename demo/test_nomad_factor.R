require(crs)

set.seed(42)

n <- 10000
num.eval <- 50

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

## No initial points, generate by sample

model.nomad <- crs(y~x1+x2,
                   basis="auto",
                   cv="nomad",
                   complexity="degree-knots",
                   knots="uniform",
                   deriv=1,
                   cv.func="cv.aic")

summary(model.nomad)

## Initial points, please note that the length of x0 should be the
## same as the number of variables in nomad. This will depend on the
## parameter complexity and the number of predictors.

x0<-c(3, 3, 2, 4)

## Single inital point

model.nomad <- crs(y~x1+x2,
                   basis="auto",
                   cv="nomad",
                   complexity="degree-knots",
                   knots="uniform",
                   deriv=1,
                   cv.func="cv.aic", 
                   x0=x0)

summary(model.nomad)

## Multiple initial points, if x0 is not provided and  best_x.txt exists,
## best_x.txt will be read in as the first inital point, otherwise, points
## will be generated randomly.

model.multi.nomad <- crs(y~x1+x2,
                         basis="auto",
                         cv="nomad",
                         complexity="degree-knots",
                         knots="uniform",
                         deriv=1,
                         cv.func="cv.aic", 
                         x0=x0, 
                         nmulti=10)

summary(model.multi.nomad)

## Compare with exhaustive search (cv="exhaustive")

model <- crs(y~x1+x2,
             basis="auto",
             cv="exhaustive",
             complexity="degree-knots",
             knots="uniform",
             deriv=1,
             cv.func="cv.aic")

summary(model)
