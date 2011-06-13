% $ID$
\name{krscvNOMAD}
\alias{krscvNOMAD}

\title{Categorical Kernel Regression Spline Cross-Validation}

\description{
  \code{krscvNOMAD} computes NOMAD-based (Nonsmooth Optimization by Mesh
Adaptive Direct Search, Abramson, Audet, Couture and Le Digabel (2011))
cross-validation directed search for a regression spline estimate of a
one (1) dimensional dependent variable on an \code{r}-dimensional vector
of continuous and ordinal/nominal factor predictors.
}
\usage{
krscvNOMAD(xz,
           y,
           degree.max = 10, 
           segments.max = 10, 
           degree.min = 0, 
           segments.min = 1, 
           kernel.type = c("nominal","ordinal"),
           restarts = 0,
           complexity = c("degree-knots","degree","knots"),
           knots = c("quantiles","uniform"),
           basis = c("additive","tensor","auto"),
           cv.func = c("cv.ls","cv.gcv","cv.aic"),
           degree = degree,
           segments = segments,
           lambda = lambda, 
           opts = list(),
           nmulti = 0)
}

\arguments{

  \item{y}{
    continuous univariate vector
  }
  
  \item{xz}{
    continuous and/or ordinal/nominal predictors
  }

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
    (\dQuote{\code{nominal}}, default) or ordered
    (\dQuote{\code{ordinal}}) - see \code{\link{crs}} for details
  }
  
\item{restarts}{
    number of times to restart search from different initial random
    values, defaults to \code{0}
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
    cross-validation. Defaults to \code{cv.ls}.  }
  
  \item{degree}{ integer/vector specifying the degree of the B-spline
  basis for each dimension of the continuous \code{x}}
  
  \item{segments}{ integer/vector specifying the number of segments of
  the B-spline basis for each dimension of the continuous \code{x}
  (i.e. number of knots minus one)}

  \item{lambda}{ real/vector for the categorical predictors. If it is not
  NULL,  it will be a starting value of lambda}

\item{opts}{
  arguments passed to the NOMAD solver.
  }
  
  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points. Defaults to \code{0}.
  }

}

\details{

  \code{krscvNOMAD} computes cross-validation for a regression spline
  estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous and ordinal/nominal factor
  predictors. The optimal \code{K}/\code{lambda} combination is returned
  along with other results (see below for return values). The method
  uses kernel functions appropriate for categorical (ordinal/nominal)
  predictors which avoids the loss in efficiency associated with
  sample-splitting procedures that are typically used when faced with a
  mix of continuous and ordinal/nominal factor predictors.

  For the continuous predictors the regression spline model employs
  either the additive or tensor product B-spline basis matrix for a
  multivariate polynomial spline via the B-spline routines in the GNU
  Scientific Library and the \code{\link{tensor.prod.model.matrix}}
  function in the \pkg{mgcv} package.

  For the discrete predictors the product kernel function is of the
  \sQuote{Li-Racine} type (see Li and Racine (2007) for details).
  
  Numerical search is undertaken using \code{\link{optim}} and the
  box-constrained \code{L-BFGS-B} method (see \code{\link{optim}} for
  details). The user may restart the algorithm as many times as desired
  via the \code{restarts} argument. The approach ascends from \code{K=0}
  through \code{degree.max} and for each value of \code{K} searches for the
  optimal bandwidths for this value of \code{K}. After the most complex
  model has been searched then the optimal \code{K}/\code{lambda}
  combination is selected. If any element of the optimal \code{K} vector
  coincides with \code{degree.max} a warning is produced and the user ought
  to restart their search with a larger value of \code{degree.max}.

  Note that, since the approach is built on top of the regression spline
  for each unique \code{z} combination, the splines must be estimable
  for all every unique \code{z} combination. When this is not the case
  the user must combine data in a manner such that this is achieved.
  
}

\value{

  \code{krscvNOMAD} returns a \code{crscv} object. Furthermore, the
    function \code{\link{summary}} supports objects of this type. The
    returned objects have the following components:

  \item{K}{ scalar/vector containing optimal degree(s) of spline or
  number of segments }
  \item{K.mat}{ vector/matrix of values of \code{K} evaluated during search }  
  \item{degree.max}{ the maximum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.max=10})}
  \item{segments.max}{ the maximum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.max=10})}
  \item{degree.min}{ the minimum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.min=0})}
  \item{segments.min}{ the minimum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.min=1})}


  \item{restarts}{ number of restarts during search, if any }
  \item{lambda}{ optimal bandwidths for categorical predictors }
  \item{lambda.mat}{ vector/matrix of optimal bandwidths for each degree of spline }
  \item{cv.func}{ objective function value at optimum }
  \item{cv.func.vec}{ vector of objective function values at each degree
  of spline or number of segments in \code{K.mat}}

}

\references{

  M.A. Abramson, C. Audet, G. Couture, J.E. Dennis, Jr., and S. Le
  Digabel (2011), \dQuote{The NOMAD project}. Software available at
  http://www.gerad.ca/nomad.

  Craven, P. and G. Wahba (1979), \dQuote{Smoothing noisy data with
  spline functions,} Numerische Mathematik, 13, 377-403.

  Hurvich, C.M. and J.S. Simonoff and C.L. Tsai (1998),
  \dQuote{Smoothing parameter selection in nonparametric regression
    using an improved Akaike information criterion,} Journal of the
  Royal Statistical Society B, 60, 271-293.

  S. Le Digabel (2011), \dQuote{Algorithm 909: NOMAD: Nonlinear
  optimization with the MADS algorithm}. ACM Transactions on
  Mathematical Software, 37(4):44:1-44:15.
  
  Ma, S. and J.S. Racine and L. Yang (2011), \dQuote{Spline Regression
  in the Presence of Categorical Predictors,} manuscript.

  Ma, S. and J.S. Racine (2011), \dQuote{Additive Regression Splines
  With Irrelevant Categorical and Continuous Regressors,} manuscript.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
    Theory and Practice,} Princeton University Press.

}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca},  Zhenghua Nie <niez@mcmaster.ca>
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

## Compute the optimal K and lambda, determine optimal number of knots, set
## spline degree for x to 3

cv <- krscvNOMAD(x=xdata,y=y,complexity="knots",degree=c(3),segments=c(5))
summary(cv)
}
\keyword{nonparametric}
