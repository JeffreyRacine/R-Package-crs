## This demo considers a setting with one endogenous regressor, one
## exogenous regressor, and one instrument. Here we might wish to use
## the plot or predict commands (see ?crsiv for details) so need to
## make sure we copy over certain components from the non-IV model.

require(crs)

## This illustration was made possible by Samuele Centorrino
## <samuele.centorrino@univ-tlse1.fr>

set.seed(42)
n <- as.numeric(readline(prompt="Input the number of observations desired: "))
method <- as.numeric(readline(prompt="Input the method (0=Landweber-Fridman, 1=Tikhonov): "))
method <- ifelse(method==0,"Landweber-Fridman","Tikhonov")
nmulti <- as.numeric(readline(prompt="Input the number of multistarts desired (e.g. 10): "))

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

y <- phi(z) + 0.2*x + u

## Sort on z (for plotting)

ivdata <- data.frame(y,z,w)
ivdata <- ivdata[order(ivdata$z),]
rm(y,z,w)
attach(ivdata)

## Note that, for plotting purposes, we need to control the value of x
## (hold non-axis variables constant). We set the value of x to its
## mean for evaluation purposes (though naturally the sample x are
## used for estimation).

model.iv <- crsiv(y=y,z=z,w=w,x=x,nmulti=nmulti,method="Landweber-Fridman")
model.noniv <- crs(y~z+x,nmulti=nmulti)

## See ?crsiv for details, but this must be done in order for crsiv
## objects to play nicely with plot and predict

model.iv$terms <- model.noniv$terms
model.iv$call <- model.noniv$call
model.iv$formula <- model.noniv$formula
model.iv$xnames <- model.noniv$xnames
model.iv$xz <- model.noniv$xz

## But now we can plot and predict with abandon...

plot(model.iv,mean=T)
