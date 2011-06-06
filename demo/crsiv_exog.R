## This demo considers a setting with one endogenous regressor, one
## exogenous regressor, and one instrument.

rm(list=ls())
library(crs)

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
## used for estimation)

model.iv <- crsiv(y=y,z=z,w=w,x=x,xeval=rep(median(x),length(x)),nmulti=nmulti,method="Landweber-Fridman")
phihat.iv <- model.iv$phihat

## Now the non-iv regression spline estimator of E(y|z), again
## controlling for the evaluation value of x.

crs.mean <- predict(crs(y~z+x,nmulti=nmulti),newdata=data.frame(z,x=rep(median(x),length(x))))

## For the plots, restrict focal attention to the bulk of the data
## (i.e. for the plotting area trim out 1/4 of one percent from each
## tail of y and z)

trim <- 0.0025

if(method=="Tikhonov")  {

  subtext <- paste("Tikhonov alpha = ",
                   formatC(model.iv$alpha,digits=3,format="fg"),
                   ", n = ", n, sep="")

} else {

  subtext <- paste("Landweber-Fridman iterations = ",
                   model.iv$num.iterations,
                   ", n = ", n,sep="")

}

curve(phi,min(z),max(z),
      xlim=quantile(z,c(trim,1-trim)),
      ylim=quantile(y,c(trim,1-trim)),
      ylab="Y",
      xlab="Z",
      main="Nonparametric Instrumental Spline Regression",
      sub=subtext,
      lwd=1,lty=1)

points(z,y,type="p",cex=.25,col="grey")

lines(z,eyz(z),lwd=1,lty=1)

lines(z,phihat.iv,col="blue",lwd=2,lty=2)

lines(z,crs.mean,col="red",lwd=2,lty=4)


legend(quantile(z,trim),quantile(y,1-trim),
       c(expression(paste(varphi(z,x),", E(y|z, x)",sep="")),
         expression(paste("Nonparametric ",hat(varphi)(z,x))),
         "Nonparametric E(y|z, x)"),
       lty=c(1,2,4),
       col=c("black","blue","red"),
       lwd=c(1,2,2))

detach(ivdata)
