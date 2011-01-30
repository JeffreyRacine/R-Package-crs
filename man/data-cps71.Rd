% $ID$
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
  Cambridge	University Press.
  
}

\keyword{datasets}

\examples{
## Example - compare local linear  smooth and regression spline for
## cps71 data

data(cps71)
attach(cps71)
require(np)

model.crs <- crs(logwage~age,complexity="degree-knots")

model.np <- npreg(logwage~age,regtype="ll")

plot(age,logwage,cex=0.25,col="grey")
lines(age,fitted(model.crs))
lines(age,fitted(model.np),col="red",lty=2)
legend(22.5,15,c("Regression Spline","Local Linear"),lty=c(1,2),col=c("black","red"))

summary(model.crs)
summary(model.np)
detach("package:np")
}
