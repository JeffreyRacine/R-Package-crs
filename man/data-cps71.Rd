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
## Example - fit a spline model for log wages as a function of age.

data(cps71)
attach(cps71)

model.crs <- crs(logwage~age,complexity="degree-knots")
plot(age,logwage,cex=0.25,col="grey",
     sub=paste("crs-CV = ", formatC(model.crs$cv.score,format="f",digits=3)))
lines(age,fitted(model.crs),lty=1,col=1)

crs.txt <- paste("crs (R-squared = ",formatC(model.crs$r.squared,format="f",digits=3),")",sep="")
legend(22.5,15,crs.txt,lty=1,col=1,bty="n")

summary(model.crs)
}
