\name{wage1}
\docType{data}
\alias{wage1}
\title{ Cross-Sectional Data on Wages }
\description{ Cross-section wage data consisting of a random sample
taken from the U.S. Current Population Survey for the year 1976. There
are 526 observations in total. }

\usage{data("wage1")}

\format{ A data frame with 24 columns, and 526 rows.

\describe{
    \item{wage}{column 1, of type \code{numeric}, average hourly earnings}
    \item{educ}{column 2, of type \code{numeric}, years of education}
    \item{exper}{column 3, of type \code{numeric}, years potential experience}
    \item{tenure}{column 4, of type \code{numeric}, years with current employer}
    \item{nonwhite}{column 5, of type \code{factor}, =\dQuote{Nonwhite} if nonwhite, \dQuote{White} otherwise}
    \item{female}{column 6, of type \code{factor}, =\dQuote{Female} if female, \dQuote{Male} otherwise}
    \item{married}{column 7, of type \code{factor}, =\dQuote{Married} if Married, \dQuote{Nonmarried} otherwise}
    \item{numdep}{column 8, of type \code{numeric}, number of dependents}
    \item{smsa}{column 9, of type \code{numeric}, =1 if live in SMSA}
    \item{northcen}{column 10, of type \code{numeric}, =1 if live in north central U.S}
    \item{south}{column 11, of type \code{numeric}, =1 if live in southern region}
    \item{west}{column 12, of type \code{numeric}, =1 if live in western region}
    \item{construc}{column 13, of type \code{numeric}, =1 if work in construc. indus.}
    \item{ndurman}{column 14, of type \code{numeric}, =1 if in nondur. manuf. indus.}
    \item{trcommpu}{column 15, of type \code{numeric}, =1 if in trans, commun, pub ut}
    \item{trade}{column 16, of type \code{numeric}, =1 if in wholesale or retail}
    \item{services}{column 17, of type \code{numeric}, =1 if in services indus.}
    \item{profserv}{column 18, of type \code{numeric}, =1 if in prof. serv. indus.}
    \item{profocc}{column 19, of type \code{numeric}, =1 if in profess. occupation}
    \item{clerocc}{column 20, of type \code{numeric}, =1 if in clerical occupation}
    \item{servocc}{column 21, of type \code{numeric}, =1 if in service occupation}
    \item{lwage}{column 22, of type \code{numeric}, log(wage)}
    \item{expersq}{column 23, of type \code{numeric}, exper\eqn{^2}{**2}}
    \item{tenursq}{column 24, of type \code{numeric}, tenure\eqn{^2}{**2}}
  }
}
\source{ Jeffrey M. Wooldridge }

\references{

  Wooldridge, J.M. (2000), \emph{Introductory Econometrics: A Modern
  Approach}, South-Western College Publishing.

}

\examples{
\dontrun{
data(wage1)

## Cross-validated model selection for spline degree and bandwidths Note
## - we override the default degree.max and segments.max here given the 
## large number of predictors and small sample size (the tensor spline basis 
## will become singular hampering search)

model <- crs(lwage~married+
             female+
             nonwhite+                
             educ+
             exper+
             tenure,
             degree.max=5, 
             segments.max=5,              
             data=wage1)
             
summary(model)

## Residual plots
plot(model)
## Partial mean plots (control for non axis predictors)
plot(model,mean=TRUE)
## Partial first derivative plots (control for non axis predictors)
plot(model,deriv=1)
## Partial second derivative plots (control for non axis predictors)
plot(model,deriv=2)
## Compare with local linear kernel regression
require(np)
model <- npreg(lwage~married+
               female+
               nonwhite+                
               educ+
               exper+
               tenure,
               regtype="ll",
               bwmethod="cv.aic",
               data=wage1)

summary(model)

## Partial mean plots (control for non axis predictors)
plot(model,common.scale=FALSE)
## Partial first derivative plots (control for non axis predictors)
plot(model,gradients=TRUE,common.scale=FALSE)
detach("package:np")
} % end dontrun
}

\keyword{datasets}
