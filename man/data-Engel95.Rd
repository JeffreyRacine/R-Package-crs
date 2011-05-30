% $ID$
\name{Engel95}
\docType{data}
\alias{Engel95}
\title{ British Family Expenditure Survey  }
\description{
British cross-sedion data consisting of a random sample taken from the British
Family Expenditure Survey. The households consist of married
couples with an employed head-of-household between the ages of 25
and 55 years. There are 1655 household-level observations in total.
}
\usage{data("Engel95")}
\format{ A data frame with 10 columns, and 1655 rows.
  \describe{

"food"     "catering" "alcohol"  "fuel"     "motor"    "fares"   
 [7] "leisure"  "logexp"   "logwages" "nkids"   
    
 \item{food}{ expenditure share on food, of type \code{numeric}}
 \item{catering}{ expenditure share on catering, of type \code{numeric}}
 \item{alcohol}{ expenditure share on alcohol, of type \code{numeric}}
 \item{fuel}{ expenditure share on fuel, of type \code{numeric}}
 \item{motor}{ expenditure share on motor, of type \code{numeric}}
 \item{fares}{ expenditure share on fares, of type \code{numeric}}
 \item{leisure}{ expenditure share on leisure, of type \code{numeric}}
 \item{logexp}{ logarithm of total expenditure, of type \code{numeric}}
 \item{logwages}{ logarithm of total earnings, of type \code{numeric}}  
  }
}
\source{ Richard Blundell and Dennis Kristensen }
\references{

  Blundell, R. and X. Chen and D. Kristensen (2007),
  \dQuote{Semi-Nonparametric IV Estimation of Shape-Invariant Engel
    Curves,} Econometrica, 75, 1613–1669.

    
  Darolles, S. and Y. Fan and J.P. Florens and E. Renault (forthcoming),
  \dQuote{Nonparametric Instrumental Regression,} Econometrica.

  Feve, F. and J.P. Florens (2010), \dQuote{The practice of
  non-parametric estimation by solving inverse problems: the example of
  transformation models,} Econometrics Journal, 13, S1-S27.


  Horowitz, J.L. (2011), \dQuote{Applied Nonparametric Instrumental
  Variables Estimation,}, Econometrica, 79, 347-394.
  
  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  
}

\keyword{datasets}

\examples{
## Example - conduct IV estimation using Landeweber Fridman iteration
## Engel95 data

data(Engel95)
## Sort on logexp (the endogenous regressor) for plotting purposes
Engel95 <- Engel95[order(Engel95$logexp),] 
attach(Engel95)

y <- food
z <- logexp
w <- logwages

k <- 10
con <- 0.5
nmulti=5

phi.mat <- matrix(NA,nrow=length(y),ncol=k)

phi.0 <- crs(y~z,nmulti=nmulti)
phi.mat[,1] <- fitted(phi.0) + fitted(crs(fitted(crs(residuals(phi.0)~w,nmulti=nmulti))~z,nmulti=nmulti))

for(j in 2:k) {
    phi.mat[,j] <- phi.mat[,j-1] + con*fitted(crs(fitted(crs((y-phi.mat[,j-1])~w,nmulti=nmulti))~z,nmulti=nmulti))
}

trim <- 0.0025

plot(z,y,
     xlim=quantile(z,c(trim,1-trim)),
     ylim=quantile(y,c(trim,1-trim)),
      ylab="Y",
      xlab="Z",
      main="Iterated Nonparametric Instrumental Variables",
     type="p",
     cex=.5,
     col="lightgrey")

lines(z,phi.mat[,k],col="blue",lwd=2,lty=2)

lines(z,fitted(phi.0),col="red",lwd=2,lty=4)

legend(quantile(z,trim),quantile(y,1-trim),
       c(expression(paste("Nonparametric IV: ",hat(varphi)(z))),
         "Nonparametric Regression: E(y|z)"),
       lty=c(2,4),
       col=c("blue","red"),
       lwd=c(2,2))


}
