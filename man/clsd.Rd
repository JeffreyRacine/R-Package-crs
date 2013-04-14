\name{clsd}
\alias{clsd}
\title{Categorical Regression Splines}

\description{ \code{clsd} computes a logspline density estimate of a one
    (1) dimensional continuous variable.}

\usage{
clsd(x=NULL,
     beta=NULL,
     xeval=NULL,
     degree=NULL,
     segments=NULL,
     degree.min=2,
     degree.max=10,
     segments.min=1,
     segments.max=10,
     basis="tensor",
     knots="quantiles",
     penalty=c("aic","sic","cv","none"),
     deriv.index=1,
     deriv=0,
     do.break=FALSE,
     er=NULL,
     linearize=TRUE,
     n.integrate=1.0e+03,
     nmulti=1)
}
\arguments{

  \item{x}{ a numeric vector of training data. }

  \item{x}{ a numeric vector of evaluation data. }  

  \item{degree}{ integer/vector specifying the polynomial degree of the
    B-spline basis for each dimension of the continuous \code{x} (default
    \code{degree=3}, i.e. cubic spline)}

  \item{segments}{ integer/vector specifying the number of segments of the
    B-spline basis for each dimension of the continuous \code{x}
    (i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier
    curve)}

  \item{degree.max}{ the maximum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.max=10})}
  
  \item{segments.max}{ the maximum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.max=10})}

  \item{degree.min}{ the minimum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.min=0})}

  \item{segments.min}{ the minimum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.min=1})}

  \item{cv.df.min}{ the minimum degrees of freedom to allow when
    conducting NOMAD-based cross-validation (default
    \code{cv.df.min=1})}

\item{basis}{ a character string (default \code{basis="auto"})
    indicating whether the additive or tensor product B-spline basis
    matrix for a multivariate polynomial spline or generalized B-spline
    polynomial basis should be used. Note this can be automatically
    determined by cross-validation if \code{cv="nomad"} or
    \code{cv="exhaustive"} and \code{basis="auto"}, and is an
    \sQuote{all or none} proposition (i.e. interaction terms for all
    predictors or for no predictors given the nature of \sQuote{tensor
    products}). Note also that if there is only one predictor this
    defaults to \code{basis="additive"} to avoid unnecessary computation
    as the spline bases are equivalent in this case }

   \item{knots}{ a character string (default \code{knots="quantiles"})
  specifying where knots are to be placed. \sQuote{quantiles} specifies
  knots placed at equally spaced quantiles (equal number of observations
  lie in each segment) and \sQuote{uniform} specifies knots placed at
  equally spaced intervals. If \code{knots="auto"}, the knot type will
  be automatically determined by cross-validation }

  \item{deriv}{ an integer \code{l} (default \code{deriv=0}) specifying
  whether to compute the univariate \code{l}th partial derivative for
  each continuous predictor (and difference in levels for each
  categorical predictor) or not and if so what order. Note that if
  \code{deriv} is higher than the spline degree of the associated
  continuous predictor then the derivative will be zero and a warning
  issued to this effect }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points (default \code{nmulti=5})
  }

  \item{penalty}{ the penalty function to use when conducting maximum
    likelihood}

  \item{do.break}{ a logical value indicating whether or not to stop
     search for a global minimum once the objective function stops
     improving}

   \item{er}{
     a scalar indicating the fraction of data range to extend
     the tails (default \code{1/log(n)})
   }

   \item{linearize}{
     a logical value indicating whether to linearize the first and last
     bases (default \code{TRUE}
   }

   \item{n.integrate}{
     the number of integration points in addition to the sample
     realizations evenly spaced on the extended range specificed by
     \code{er} (defaults to \code{1.0e+03})
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
    \code{\link{fitted}} and \code{\link{residuals}} extract (or
    generate) estimated values and residuals. Furthermore, the functions
    \code{\link{summary}}, \code{\link{predict}}, and \code{\link{plot}}
    (options \code{mean=FALSE}, \code{deriv=i} where \eqn{i} is an
    integer, \code{ci=FALSE}, \code{persp.rgl=FALSE},
    \code{plot.behavior=c("plot","plot-data","data")},
    \code{xtrim=0.0},\code{xq=0.5}) support objects of this type. The
    returned object has the following components:

  \item{fitted.values}{ estimates of the regression function
    (conditional mean) at the sample points or evaluation points }

  \item{degree}{ integer/vector specifying the degree of the B-spline
    basis for each dimension of the continuous \code{x}}
  
  \item{segments}{ integer/vector specifying the number of segments of
  the B-spline basis for each dimension of the continuous \code{x}}

  \item{deriv.mat}{ a matrix of derivatives (or differences in levels
  for the categorical \code{z}) whose order is determined by
  \code{deriv=} in the \code{clsd} call }

}
  
\section{Usage Issues}{


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
set.seed(42)
require(logspline)
## Example - simulated data
n <- 1000
df <- 2
x <- sort(rchisq(n,df=df))
f.dgp <- dchisq(x,df=df)

clsd.out <- clsd(x)
f.hat <- clsd.out$density
## Standard (cubic) estimate taken from the logspline package
f.lsc <- dlogspline(x,fit.cub<-logspline(x))
## Compute MSEs
mse.f.hat <- mean((f.hat-f.dgp)^2)
mse.f.lsc <- mean((f.lsc-f.dgp)^2)
## Compute ylimits for plot
ylim <- c(0,max(f.hat,f.lsc,f.dgp))

plot(x,f.dgp,
     ylim=ylim,
     type="l",
     ylab="Density",
     xlab="X",
     sub=paste("MSE: kernel = ",format(mse.f.kernel),", Cubic = ",format(mse.f.lsc),", Proposed = ", format(mse.f.hat)))

lines(x,f.lsc,col=2,lty=2)
lines(x,f.hat,col=3,lty=1)

rug(x)

legend("topright",c("DGP",paste("Cubic Logspline Density (knots = ",fit.cub$nknots,")",sep=""),
                    paste("Logspline Density (degree = ", clsd.out$degree, ", segments = ",clsd.out$segments,", cv = ",scan("penalty.dat",what="char"),")",sep="")),
       lty=c(1:3,1),
       col=1:3,
       bty="n",
       cex=0.6)
} %% End dontrun
}
\keyword{nonparametric}

