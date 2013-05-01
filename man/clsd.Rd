\name{clsd}
\alias{clsd}
\title{Categorical Logspline Density}

\description{ \code{clsd} computes the logspline density, density
    derivative, distribution, and smoothed quantiles for a one (1)
    dimensional continuous variable using the approach of Racine
    (2013).}

\usage{
clsd(x = NULL,
     beta = NULL,
     xeval = NULL,
     degree = NULL,
     segments = NULL,
     degree.min = 2,
     degree.max = 10,
     segments.min = 1,
     segments.max = 10,
     lbound = NULL,
     ubound = NULL,
     basis = "tensor",
     knots = "quantiles",
     penalty = NULL,
     deriv.index = 1,
     deriv = 1,
     elastic.max = TRUE,
     elastic.diff = 3,
     do.gradient = TRUE,
     er = NULL,
     monotone = TRUE,
     monotone.lb = -250,
     n.integrate = 500,
     nmulti = 1,
     method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
     verbose = FALSE,
     quantile.seq = seq(.01,.99,by=.01),
     random.seed = 42,
     maxit = 10^5,
     max.attempts = 25,
     NOMAD = FALSE)
}
\arguments{

  \item{x}{ a numeric vector of training data }

  \item{beta}{ a numeric vector of coefficients (default \code{NULL}) }

  \item{xeval}{ a numeric vector of evaluation data }  

  \item{degree}{ integer/vector specifying the polynomial degree of the
    B-spline basis for each dimension of the continuous \code{x} (default
    \code{degree=2})}

  \item{segments}{ integer/vector specifying the number of segments of the
    B-spline basis for each dimension of the continuous \code{x}
    (i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier
    curve)}

  \item{segments.min,segments.max}{ when \code{elastic.max=FALSE}, the
    minimum/maximum segments of the B-spline basis for each of the
    continuous predictors (default
    \code{segments.min=1},\code{segments.max=10})}

  \item{degree.min,degree.max}{ when \code{elastic.max=FALSE} the
    minimum/maximum degree of the B-spline basis for each of the
    continuous predictors (default \code{degree.min=2},
    \code{degree.max=10})}

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

  \item{deriv}{ an integer \code{l} (default \code{deriv=1}) specifying
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
    points (default \code{nmulti=1})
  }

  \item{penalty}{ the parameter to be used in the AIC criterion. The
    method chooses the number of degrees plus number of segments
    (knots-1) that maximizes \code{2*logl-penalty*(degree+segments)}. The
    default is to use the penalty parameter of \code{log(n)/2} (\code{2}
    would deliver standard AIC, \code{log(n)} standard BIC)}
  
  \item{elastic.max,elastic.diff}{ a logical value/integer indicating
     whether to use `elastic' search bounds such that the optimal
     degree/segment must lie \code{elastic.diff} units from the
     respective search bounds }

  \item{do.gradient}{ a logical value indicating whether or not to use
     the analytical gradient during optimization (defaults to \code{TRUE})}

   \item{er}{
     a scalar indicating the fraction of data range to extend
     the tails (default \code{1/log(n)}, see \code{\link{extendrange}} for
     further details)
   }

   \item{monotone}{ a logical value indicating whether modify
     the standard B-spline basis function so that it is tailored for
     density estimation (default \code{TRUE})
   }

   \item{monotone.lb}{ a negative bound specifying the lower bound on
     the linear segment coefficients used when (\code{monotone=FALSE})
   }

   \item{n.integrate}{
     the number of evenly spaced integration points on the extended range specified by \code{er} (defaults to \code{500})
   }

   \item{method}{
     see \code{\link{optim}} for details
   }

   \item{verbose}{
     a logical value which when \code{TRUE} produces verbose output
     during optimization
   }

   \item{quantile.seq}{
     a sequence of numbers lying in \eqn{[0,1]} on which quantiles from
     the logspline distribution are obtained
   }

   \item{random.seed}{ seeds the random number generator for initial
     parameter values when \code{\link{optim}} is called }

  \item{maxit}{
    maximum number of iterations used by \code{\link{optim}}
   }
   
   \item{max.attempts}{
     maximum number of attempts to undertake if \code{\link{optim}}
     fails for any set of initial parameters for each value of
     \code{nmulti}
   }

   \item{NOMAD}{
     a logical value which when \code{TRUE} calls \code{\link{snomadr}}
     to determine the optimal \code{degree} and \code{segments}
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
    \code{\link{coef}}, \code{\link{fitted}}, \code{\link{plot}} and
    \code{\link{summary}} support objects of this type (\code{er=FALSE}
    plots the density on the sample realizations (default is `extended
    range' data), see \code{er} above, \code{distribution=TRUE} plots
    the distribution). The returned object has the following components:

  \item{density}{ estimates of the density function
    at the sample points}

  \item{density.er}{ the density evaluated on the `extended range'
  of the data }

  \item{density.deriv}{ estimates of the derivative of the density function
    at the sample points}

  \item{density.deriv.er}{ estimates of the derivative of the density
    function evaluated on the `extended range' of the data }

  \item{distribution}{ estimates of the distribution function
    at the sample points}

  \item{distribution.er}{ the distribution evaluated on the `extended range'
  of the data }

  \item{xer}{ the `extended range' of the data }

  \item{degree}{ integer/vector specifying the degree of the B-spline
    basis for each dimension of the continuous \code{x}}
  
  \item{segments}{ integer/vector specifying the number of segments of
    the B-spline basis for each dimension of the continuous \code{x}}

  \item{xq}{ vector of quantiles }

  \item{tau}{ vector generated by \code{quantile.seq} or input by the
  user (lying in \code{[0,1]}) from which the quantiles \code{xq} are
  obtained}

}
  
\section{Usage Issues}{

  This function should be considered to be in `beta' status until
  further notice.

  If smoother estimates are desired and \code{degree=degree.min}, increase
  \code{degree.min} to, say, \code{degree.min=3}.

  The use of `regression' B-splines can lead to undesirable behavior at
  the endpoints of the data (i.e. when \code{monotone=FALSE}). The
  default `density' B-splines ought to be well-behaved in these regions.

}
\references{

  Racine, J.S. (2013), \dQuote{Logspline Mixed Data Density Estimation,}
  manuscript.

}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\seealso{
  \code{\link[logspline]{logspline}}
}

\examples{
\dontrun{
## Old Faithful eruptions data histogram and clsd density

library(MASS)
data(faithful)
attach(faithful)

model <- clsd(eruptions)

ylim <- c(0,max(model$density,hist(eruptions,breaks=20,plot=FALSE)$density))

plot(model,ylim=ylim)

hist(eruptions,breaks=20,freq=FALSE,add=TRUE,lty=2)

rug(eruptions)

summary(model)

coef(model)

## Simulated data

set.seed(42)
require(logspline)

## Example - simulated data

n <- 250
x <- sort(rnorm(n))
f.dgp <- dnorm(x)

model <- clsd(x)

## Standard (cubic) estimate taken from the logspline package
## Compute MSEs

mse.clsd <- mean((fitted(model)-f.dgp)^2)

model.logspline <- logspline(x)

mse.logspline <- mean((dlogspline(x,model.logspline)-f.dgp)^2)

ylim <- c(0,max(fitted(model),dlogspline(x,model.logspline),f.dgp))

plot(model,
     ylim=ylim,
     sub=paste("MSE: logspline = ",format(mse.logspline),", clsd = ",
     format(mse.clsd)),
     lty=3,
     col=3)

xer <- model$xer

lines(xer,dlogspline(xer,model.logspline),col=2,lty=2)
lines(xer,dnorm(xer),col=1,lty=1)

rug(x)

legend("topright",c("DGP",paste("Cubic Logspline Density (package `logspline', knots = ",model.logspline$nknots,")",sep=""),
                    paste("clsd Density (degree = ", model$degree, ", segments = ",model$segments,", penalty = ",round(model$penalty,2),")",sep="")),
       lty=1:3,
       col=1:3,
       bty="n",
       cex=0.75)

summary(model)

coef(model)

## Simulate data with known bounds

set.seed(42)
n <- 10000
x <- runif(n,0,1)

model <- clsd(x,lbound=0,ubound=1)

plot(model)
} %% End dontrun
}
\keyword{nonparametric}

