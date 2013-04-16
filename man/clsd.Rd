\name{clsd}
\alias{clsd}
\title{Categorical Logspline Density}

\description{ \code{clsd} computes a logspline density estimate of a one
    (1) dimensional continuous variable.}

\usage{
clsd(x=NULL,
     beta=NULL,
     xeval=NULL,
     degree=NULL,
     segments=NULL,
     degree.min=2,
     degree.max=5,
     segments.min=1,
     segments.max=5,
     lbound=NULL,
     ubound=NULL,
     basis="tensor",
     knots="quantiles",
     penalty=c("aic","sic","cv","none"),
     deriv.index=1,
     deriv=0,
     do.break=FALSE,
     do.gradient=TRUE,
     er=NULL,
     monotone=TRUE,
     n.integrate=1.0e+03,
     nmulti=1,
     method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
     verbose=FALSE)
}
\arguments{

  \item{x}{ a numeric vector of training data }

  \item{beta}{ a numeric vector of coefficients (default \code{NULL}) }

  \item{xeval}{ a numeric vector of evaluation data }  

  \item{degree}{ integer/vector specifying the polynomial degree of the
    B-spline basis for each dimension of the continuous \code{x} (default
    \code{degree=3}, i.e. cubic spline)}

  \item{segments}{ integer/vector specifying the number of segments of the
    B-spline basis for each dimension of the continuous \code{x}
    (i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier
    curve)}

  \item{segments.min,segments.max}{ the minimum/maximum segments of the
    B-spline basis for each of the continuous predictors (default
    \code{segments.min=1},\code{segments.max=5})}

  \item{degree.min,degree.max}{ the minimum/maximum degree of the
    B-spline basis for each of the continuous predictors (default
    \code{degree.min=2}, \code{degree.max=5})}

  \item{lbound,ubound}{
    lower/upper bound for the support of the density. For example, if
    there is a priori knowledge that the density equals zero to the left
    of 0, and has a discontinuity at 0, the user could specify lbound =
    0. However, if the density is essentially zero near 0, one does not
    need to specify lbound }

  \item{basis}{ a character string (default \code{basis="tensor"})
    indicating whether the additive or tensor product B-spline basis
    matrix for a multivariate polynomial spline or generalized B-spline
    polynomial basis should be used }

   \item{knots}{ a character string (default \code{knots="quantiles"})
  specifying where knots are to be placed. \sQuote{quantiles} specifies
  knots placed at equally spaced quantiles (equal number of observations
  lie in each segment) and \sQuote{uniform} specifies knots placed at
  equally spaced intervals }

  \item{deriv}{ an integer \code{l} (default \code{deriv=0}) specifying
  whether to compute the univariate \code{l}th partial derivative for
  each continuous predictor (and difference in levels for each
  categorical predictor) or not and if so what order. Note that if
  \code{deriv} is higher than the spline degree of the associated
  continuous predictor then the derivative will be zero and a warning
  issued to this effect }

  \item{deriv.index}{ an integer \code{l} (default \code{deriv.index=1})
  specifying the index (currently only supports 1) of the variable whose
  derivative is requested  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points (default \code{nmulti=5})
  }

  \item{penalty}{ the penalty function to use when conducting maximum
    likelihood}

  \item{do.break}{ a logical value indicating whether or not to stop
     search for a global maximum once the objective function stops
     improving}

  \item{do.gradient}{ a logical value indicating whether or not to use
     the analytical gradient during optimization }

   \item{er}{
     a scalar indicating the fraction of data range to extend
     the tails (default \code{1/log(n)}, see \link{extendrange} for
     further details)
   }

   \item{monotone}{ a logical value indicating whether modify
     the standard B-spline basis function so that it is better suited to
     density estimation (default \code{TRUE})
   }

   \item{n.integrate}{
     the number of integration points in addition to the sample
     realizations evenly spaced on the extended range specified by
     \code{er} (defaults to \code{1.0e+03})
   }

   \item{method}{
     see \link{optim} for details
   }

   \item{verbose}{
     a logical value which when \code{TRUE} produces verbose output
     during optimization
   }

  
}

\details{

  Typical usages are (see below for a  list of options and also
  the examples at the end of this help file)
  \preformatted{
    
    model <- clsd(x)

    }
    
  \code{clsd} computes a logspline density estimate of a one (1)
  dimensional continuous variable.

  The spline model employs the tensor product B-spline basis matrix for
  a multivariate polynomial spline via the B-spline routines in the GNU
  Scientific Library (\url{http://www.gnu.org/software/gsl/}) and the
  \code{\link{tensor.prod.model.matrix}} function.

  When \code{basis="additive"} the model becomes additive in nature
  (i.e. no interaction/tensor terms thus semiparametric not fully
  nonparametric).

  When \code{basis="tensor"} the model uses the multivariate tensor
  product basis.

}

\value{

  \code{clsd} returns a \code{clsd} object.  The generic functions
    \code{\link{summary}} and \code{\link{plot}} support objects of this
    type (\code{er=TRUE} plots the density on `extended range' data, see
    \code{er} above, \code{distribution=TRUE} plots the
    distribution). The returned object has the following components:

  \item{density}{ estimates of the density function
    at the sample points}

  \item{density.xer}{ the density evaluated on the `extended range'
  of the data }

  \item{density.deriv}{ estimates of the derivative of the density function
    at the sample points}

  \item{distribution}{ estimates of the distribution function
    at the sample points}

  \item{distribution.xer}{ the distribution evaluated on the `extended range'
  of the data }

  \item{xer}{ the `extended range' of the data }

  \item{degree}{ integer/vector specifying the degree of the B-spline
    basis for each dimension of the continuous \code{x}}
  
  \item{segments}{ integer/vector specifying the number of segments of
  the B-spline basis for each dimension of the continuous \code{x}}

}
  
\section{Usage Issues}{

  The use of a single linear segment for the left and rightmost spline
  bases can have lead to undesirable tail behavior beyond the range of
  the data when \code{monotone=FALSE} (the standard B-spline bases are
  vertical at the min/max of \code{x} which can lead to kinks in the
  estimated density at these points).

}
\references{

  Ma, S. and J.S. Racine and L. Yang (2011), \dQuote{Spline Regression
  in the Presence of Categorical Predictors,} manuscript.

  Ma, S. and J.S. Racine (forthcoming), \dQuote{Additive Regression
  Splines with Irrelevant Categorical and Continuous Regressors,}
  Statistica Sinica.

}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

%\section{Usage Issues}{
%}

\seealso{
  \code{\link[logspline]{logspline}}
}

\examples{
\dontrun{
## Simulated data
set.seed(42)
require(logspline)
## Example - simulated data
n <- 250
x <- sort(rnorm(n))
f.dgp <- dnorm(x)

clsd.out <- clsd(x,degree.max=3,segments.max=3)
f.clsd <- clsd.out$density
## Standard (cubic) estimate taken from the logspline package
f.lsc <- dlogspline(x,fit.cub<-logspline(x))
## Compute MSEs
mse.f.clsd <- mean((f.clsd-f.dgp)^2)
mse.f.lsc <- mean((f.lsc-f.dgp)^2)
## Compute ylimits for plot
ylim <- c(0,max(f.clsd,f.lsc,f.dgp))

plot(x,f.dgp,
     ylim=ylim,
     type="l",
     ylab="Density",
     xlab="X",
     sub=paste("MSE: logspline = ",format(mse.f.lsc),", clsd = ", format(mse.f.clsd)))

lines(x,f.lsc,col=2,lty=2)
lines(x,f.clsd,col=3,lty=3)

rug(x)

legend("topright",c("DGP",paste("Cubic Logspline Density (package `logspline', knots = ",fit.cub$nknots,")",sep=""),
                    paste("clsd Density (degree = ", clsd.out$degree, ", segments = ",clsd.out$segments,", cv = ",clsd.out$penalty,")",sep="")),
       lty=1:3,
       col=1:3,
       bty="n",
       cex=0.75)

## Old Faithful eruptions data histogram and clsd density
library(MASS)
data(faithful)
attach(faithful)
clsd.out <- clsd(eruptions)
ylim <- c(0,max(clsd.out$density,hist(eruptions,breaks=20,plot=FALSE)$density))
plot(clsd.out,ylim=ylim)
hist(eruptions,breaks=20,freq=FALSE,add=TRUE,lty=2)
rug(eruptions)
summary(clsd.out)
} %% End dontrun
}
\keyword{nonparametric}

