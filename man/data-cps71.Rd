\name{cps71}
\docType{data}
\alias{cps71}
\title{ Canadian High School Graduate Earnings  }
\description{
  Canadian cross-section wage data consisting of a random sample taken
  from the 1971 Canadian Census Public Use Tapes for male individuals
  having common education (grade 13). There are 205 observations in total.
}
\usage{data("cps71")}
\format{ A data frame with 2 columns, and 205 rows.
  \describe{
    \item{logwage}{ the first column, of type \code{numeric}}
    \item{age}{ the second column, of type \code{integer}}
  }
}
\source{ Aman Ullah }
\references{
  
  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.
  
}

\keyword{datasets}

\examples{
## Example - we compare the nonparametric local linear kernel regression
## method with the regression spline for the cps71 data. Note that there
## are no categorical predictors in this dataset so we are merely
## comparing and contrasting the two nonparametric estimates.

data(cps71)
attach(cps71)
require(np)

model.crs <- crs(logwage~age,complexity="degree-knots")
model.np <- npreg(logwage~age,regtype="ll")

plot(age,logwage,cex=0.25,col="grey",
     sub=paste("crs-CV = ", formatC(model.crs$cv.score,format="f",digits=3),
       ", npreg-CV = ", formatC(model.np$bws$fval,format="f",digits=3),sep=""))
lines(age,fitted(model.crs),lty=1,col=1)
lines(age,fitted(model.np),lty=2,col=2)

crs.txt <- paste("crs (R-squared = ",formatC(model.crs$r.squared,format="f",digits=3),")",sep="")
np.txt <- paste("ll-npreg (R-squared = ",formatC(model.np$R2,format="f",digits=3),")",sep="")

legend(22.5,15,c(crs.txt,np.txt),lty=c(1,2),col=c(1,2))

summary(model.crs)
summary(model.np)
detach("package:np")
}
