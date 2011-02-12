\name{gsl.bs}
\alias{gsl.bs}
\alias{gsl.bs.default}
\alias{gsl.bs.predict}
\title{GSL (GNU Scientific Library) B-spline/B-spline Derivatives}

\description{ \code{gsl.bs} generates the B-spline basis matrix for a
  polynomial spline }

\usage{
gsl.bs(\dots)
\method{gsl.bs}{default}(x,
       degree = 3,
       nbreak = 2,
       deriv = 0,
       x.min = NULL,
       x.max = NULL,
       intercept = FALSE,
       knots = NULL,
       \dots)
}
\arguments{

\item{x}{ the predictor variable.  Missing values are not allowed }

\item{degree}{ degree of the piecewise polynomial - default is
\sQuote{3} (cubic spline) }
        
\item{nbreak}{ number of breaks in each interval - default is \sQuote{2}}

\item{deriv}{ the order of the derivative to be computed-default if
\code{0}}

\item{x.min}{ the lower bound on which to construct the spline -
defaults to \code{min(x)}}
        
\item{x.max}{ the upper bound on which to construct the spline -
defaults to \code{max(x)}}
        
\item{intercept}{ if \sQuote{TRUE}, an intercept is included in the
basis; default is \sQuote{FALSE} }

\item{knots}{ a vector (default \code{knots="NULL"}) specifying knots
  for the spline basis (default enables uniform knots, otherwise those
  provided are used)}

\item{\dots}{ optional arguments }
  
}

\details{

  Typical usages are (see below for a  list of options and also
  the examples at the end of this help file)
  \preformatted{
    B <- gsl.bs(x,degree=10)
    }
    
}

\value{

  \code{gsl.bs} returns a \code{gsl.bs} object.  A matrix of dimension
  \sQuote{c(length(x), degree)}.

  A primary use is in modelling formulas to directly specify a piecewise
     polynomial term in a model.  }
  
\references{

  Ma, S. and J.S. Racine and L. Yang (2010), \dQuote{Spline Regression
  in the Presence of Categorical Predictors,} manuscript.

  Ma, S. and J.S. Racine (2011), \dQuote{Additive Regression Splines
  With Irrelevant Categorical and Continuous Regressors,} manuscript.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
		Theory and Practice,} Princeton University Press.

}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

%\section{Usage Issues}{
%}

\seealso{
  \code{\link[splines]{bs}}
}

\examples{
## Plot the spline bases and their first order derivatives
x <- seq(0,1,length=100)
matplot(x,gsl.bs(x,degree=5),type="l")
matplot(x,gsl.bs(x,degree=5,deriv=1),type="l")

## Regression example
n <- 1000
x <- sort(runif(n))
y <- cos(2*pi*x) + rnorm(n,sd=.25)
B <- gsl.bs(x,degree=5,intercept=FALSE)
plot(x,y,cex=.5,col="grey")
lines(x,fitted(lm(y~B)))
}
\keyword{nonparametric}

