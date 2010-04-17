% $ID$
\name{krscv}
\alias{krscv}

\title{Categorical Kernel Regression Spline Cross-Validation}

\description{
  
  \code{krscv} computes cross-validation for regression spline
  estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous and ordinal/nominal factor
  predictors.
  
}

\usage{
krscv(xz,
      y,
      max.K=10,
      kernel.type=c("nominal","ordinal"),
      restarts=0,
      tensor = c("enabled","disabled","auto"),
      cv.norm = c("L2","L1"))
}

\arguments{

  \item{y}{
    continuous univariate vector
  }
  
  \item{xz}{
    continuous and/or ordinal/nominal predictors
  }

  \item{max.K}{
    integer specifying the maximum degree of the regression spline,
    defaults to \code{10}
  }

  \item{kernel.type}{
    a character string indicating whether the kernel
    for the discrete predictors \code{z} should be unordered
    (\dQuote{\code{nominal}}, default) or ordered
    (\dQuote{\code{ordinal}}) - see \code{\link{crs}} for details
  }
  
\item{restarts}{
    number of times to restart search from different initial random
    values, defaults to \code{0}
  }

  \item{tensor}{ a character string (default \code{tensor="enabled"})
  indicating whether tensor product terms (i.e. interaction terms) are
  included in the model or not. Note this is automatically determined
  when \code{tensor="auto"}, and is an \sQuote{all or none}
  proposition (i.e. interaction terms for all predictors or for no
  predictors) }

  \item{cv.norm}{ a character string (default \code{cv.norm="L2"}
  indicating whether the L2-norm is used (sum of squared delete-one
  residuals) or the L1-norm (sum of absolute delete-one-residuals) in
  the cross-validation function }
}

\details{

  \code{krscv} computes cross-validation for a regression spline
  estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous and ordinal/nominal factor
  predictors. The optimal \code{K}/\code{lambda} combination is returned
  along with other results (see below for return values). The method
  uses kernel functions appropriate for categorical (ordinal/nominal)
  predictors which avoids the loss in efficiency associated with
  sample-splitting procedures that are typically used when faced with a
  mix of continuous and ordinal/nominal factor predictors.

  For the continuous predictors the regression spline model employs the
  tensor product B-spline basis matrix for a multivariate polynomial
  spline via the B-spline routines in the GNU Scientific Library and the
  \code{\link{tensor.prod.model.matrix}} function in the \pkg{mgcv}
  package.

  For the discrete predictors the product kernel function is of the
  \sQuote{Li-Racine} type (see Li and Racine (2007) for details).
  
  Numerical search is undertaken using \code{\link{optim}} and the
  box-constrained \code{L-BFGS-B} method (see \code{\link{optim}} for
  details). The user may restart the algorithm as many times as desired
  via the \code{restarts} argument. The approach ascends from \code{K=0}
  through \code{max.K} and for each value of \code{K} searches for the
  optimal bandwidths for this value of \code{K}. After the most complex
  model has been searched then the optimal \code{K}/\code{lambda}
  combination is selected. If any element of the optimal \code{K} vector
  coincides with \code{max.K} a warning is produced and the user ought
  to restart their search with a larger value of \code{max.K}.

  Note that, since the approach is built on top of the regression spline
  for each unique \code{z} combination, the splines must be estimable
  for all every unique \code{z} combination. When this is not the case
  the user must combine data in a manner such that this is achieved.
  
}

\value{

  \code{krscv} returns a \code{crscv} object. Furthermore, the
    function \code{\link{summary}} supports objects of this type. The
    returned objects have the following components:

  \item{K}{ scalar/vector containing optimal degree(s) of spline }
  \item{K.mat}{ vector/matrix of values of \code{K} evaluated during search }  
  \item{max.K}{ maximum degree of spline }
  \item{restarts}{ number of restarts during search, if any }
  \item{lambda}{ optimal bandwidths for categorical predictors }
  \item{lambda.mat}{ vector/matrix of optimal bandwidths for each degree of spline }
  \item{cv.func}{ objective function value at optimum }
  \item{cv.func.vec}{ vector of objective function values at each degree
  of spline in \code{K.mat}}

}

\references{

  Racine, J.S and L. Yang and S. Ma (2010), \dQuote{Spline Regression in
  the Presence of Categorical Predictors,} manuscript.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
		Theory and Practice,} Princeton University Press.

}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

%\section{Usage Issues}{
%}

\seealso{
  \code{\link{loess}}, \code{\link[np]{npregbw}}, 
}

\examples{
set.seed(123)
## Simulated data
n <- 10000

x <- runif(n)
z <- round(runif(n,min=-0.5,max=1.5))
z.unique <- uniquecombs(as.matrix(z))
ind <-  attr(z.unique,"index")
ind.vals <-  sort(unique(ind))
dgp <- numeric(length=n)
for(i in 1:nrow(z.unique)) {
  zz <- ind == ind.vals[i]
  dgp[zz] <- z[zz]+cos(2*pi*x[zz])
}
y <- dgp + rnorm(n,sd=.1)

xdata <- data.frame(x,z=factor(z))

## Compute the optimal K and lambda

cv <- krscv(x=xdata,y=y)
summary(cv)
}
\keyword{nonparametric}
