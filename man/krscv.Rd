\name{krscv}
\alias{krscv}

\title{Categorical Kernel Regression Spline Cross-Validation}

\description{
  
  \code{krscv} computes exhaustive cross-validation directed search for
  a regression spline estimate of a one (1) dimensional dependent
  variable on an \code{r}-dimensional vector of continuous and
  nominal/ordinal (\code{\link{factor}}/\code{\link{ordered}})
  predictors.
  
}

\usage{
krscv(xz,
      y,
      degree.max=10, 
      segments.max=10,
      degree.min=0,
      segments.min=1, 
      kernel.type = c("nominal","ordinal"),
      restarts = 0,
      complexity = c("degree-knots","degree","knots"),
      knots = c("quantiles","uniform"),
      basis = c("additive","tensor","auto"),
      cv.func = c("cv.ls","cv.gcv","cv.aic"),
      degree = degree,
      segments = segments)
}

\arguments{

  \item{y}{
    continuous univariate vector
  }
  
  \item{xz}{ continuous and/or nominal/ordinal
    (\code{\link{factor}}/\code{\link{ordered}}) predictors }

  \item{degree.max}{ the maximum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.max=10})}
  
  \item{segments.max}{ the maximum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.max=10})}
  
  \item{degree.min}{ the minimum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.min=0})}
  
  \item{segments.min}{ the minimum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.min=1})}

  \item{kernel.type}{
    a character string indicating whether the kernel
    for the discrete predictors \code{z} should be unordered
    (default \dQuote{\code{nominal}}) or ordered
    (\dQuote{\code{ordinal}}) - see \code{\link{crs}} for details
  }
  
  \item{restarts}{
    number of times to restart \code{\link{optim}} from different initial
    random values (default \code{restarts=0}) when searching for optimal
    bandwidths for the categorical predictors for each unique \code{K}
    combination (i.e.\ \code{degree}/\code{segments})
  }

  \item{complexity}{a character string (default
  \code{complexity="degree-knots"}) indicating whether model `complexity' is
  determined by the degree of the spline or by the number of segments
  (`knots'). This option allows the user to use cross-validation to
  select either the spline degree (number of knots held fixed) or the
  number of knots (spline degree held fixed) or both the spline degree
  and number of knots}

   \item{knots}{ a character string (default \code{knots="quantiles"})
  specifying where knots are to be placed. \sQuote{quantiles} specifies
  knots placed at equally spaced quantiles (equal number of observations
  lie in each segment) and \sQuote{uniform} specifies knots placed at
  equally spaced intervals }

  \item{basis}{ a character string (default \code{basis="additive"})
    indicating whether the additive or tensor product B-spline basis
    matrix for a multivariate polynomial spline should be used. Note
    this can be automatically determined by cross-validation if
    \code{cv=TRUE} and \code{basis="auto"}, and is an \sQuote{all or
    none} proposition (i.e. interaction terms for all predictors or for
    no predictors given the nature of \sQuote{tensor products}). Note
    also that if there is only one predictor this defaults to
    \code{basis="additive"} to avoid unnecessary computation as the
    spline bases are equivalent in this case }

  \item{cv.func}{a character string (default \code{cv.func="cv.ls"})
    indicating which method to use to select smoothing
    parameters. \code{cv.gcv} specifies generalized cross-validation
    (Craven and Wahba (1979)), \code{cv.aic} specifies expected
    Kullback-Leibler cross-validation (Hurvich, Simonoff, and Tsai
    (1998)), and \code{cv.ls} specifies least-squares
    cross-validation  }
  
  \item{degree}{ integer/vector specifying the degree of the B-spline
  basis for each dimension of the continuous \code{x}}
  
  \item{segments}{ integer/vector specifying the number of segments of
  the B-spline basis for each dimension of the continuous \code{x}
  (i.e. number of knots minus one)}

}

\details{

  \code{krscv} computes exhaustive cross-validation for a regression
  spline estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous and nominal/ordinal
  (\code{\link{factor}}/\code{\link{ordered}}) predictors. The optimal
  \code{K}/\code{lambda} combination is returned along with other
  results (see below for return values). The method uses kernel
  functions appropriate for categorical (ordinal/nominal) predictors
  which avoids the loss in efficiency associated with sample-splitting
  procedures that are typically used when faced with a mix of continuous
  and nominal/ordinal (\code{\link{factor}}/\code{\link{ordered}})
  predictors.

  For the continuous predictors the regression spline model employs
  either the additive or tensor product B-spline basis matrix for a
  multivariate polynomial spline via the B-spline routines in the GNU
  Scientific Library (\url{http://www.gnu.org/software/gsl/}) and the
  \code{\link{tensor.prod.model.matrix}} function.

  For the discrete predictors the product kernel function is of the
  \sQuote{Li-Racine} type (see Li and Racine (2007) for details).
  
  For each unique combination of \code{degree} and \code{segment},
  numerical search for the bandwidth vector \code{lambda} is undertaken
  using \code{\link{optim}} and the box-constrained \code{L-BFGS-B}
  method (see \code{\link{optim}} for details). The user may restart the
  \code{\link{optim}} algorithm as many times as desired via the
  \code{restarts} argument. The approach ascends from \code{K=0} through
  \code{degree.max}/\code{segments.max} and for each value of \code{K}
  searches for the optimal bandwidths for this value of \code{K}. After
  the most complex model has been searched then the optimal
  \code{K}/\code{lambda} combination is selected. If any element of the
  optimal \code{K} vector coincides with
  \code{degree.max}/\code{segments.max} a warning is produced and the
  user ought to restart their search with a larger value of
  \code{degree.max}/\code{segments.max}.

}

\value{

  \code{krscv} returns a \code{crscv} object. Furthermore, the
    function \code{\link{summary}} supports objects of this type. The
    returned objects have the following components:

  \item{K}{ scalar/vector containing optimal degree(s) of spline or
  number of segments }
  \item{K.mat}{ vector/matrix of values of \code{K} evaluated during search }  
  \item{restarts}{ number of restarts during search, if any }
  \item{lambda}{ optimal bandwidths for categorical predictors }
  \item{lambda.mat}{ vector/matrix of optimal bandwidths for each degree of spline }
  \item{cv.func}{ objective function value at optimum }
  \item{cv.func.vec}{ vector of objective function values at each degree
  of spline or number of segments in \code{K.mat}}

}

\references{

  Craven, P. and G. Wahba (1979), \dQuote{Smoothing noisy data with
  spline functions,} Numerische Mathematik, 13, 377-403.

  Hurvich, C.M. and J.S. Simonoff and C.L. Tsai (1998),
  \dQuote{Smoothing parameter selection in nonparametric regression
    using an improved Akaike information criterion,} Journal of the
  Royal Statistical Society B, 60, 271-293.

  Ma, S. and J.S. Racine and L. Yang (2011), \dQuote{Spline Regression
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
  \code{\link{loess}}, \code{\link[np]{npregbw}}, 
}

\examples{
set.seed(42)
## Simulated data
n <- 1000

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

## Compute the optimal K and lambda, determine optimal number of knots, set
## spline degree for x to 3

cv <- krscv(x=xdata,y=y,complexity="knots",degree=c(3))
summary(cv)
}
\keyword{nonparametric}
