\name{krscvNOMAD}
\alias{krscvNOMAD}

\title{Categorical Kernel Regression Spline Cross-Validation}

\description{
  
  \code{krscvNOMAD} computes NOMAD-based (Nonsmooth Optimization by Mesh
  Adaptive Direct Search, Abramson, Audet, Couture and Le Digabel
  (2011)) cross-validation directed search for a regression spline
  estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous and nominal/ordinal
  (\code{\link{factor}}/\code{\link{ordered}}) predictors.
  
}
\usage{
krscvNOMAD(xz,
           y,
           basis = c("additive","tensor","glp","auto"),
           complexity = c("degree-knots","degree","knots"),
           cv.df.min = 1,
           cv.func = c("cv.ls","cv.gcv","cv.aic"),
           degree = degree,
           degree.max = 10, 
           degree.min = 0, 
           display.nomad.progress = TRUE,
           display.warnings = TRUE,
           initial.mesh.size.integer = "1",
           initial.mesh.size.real = "r0.1",
           knots = c("quantiles","uniform","auto"),
           lambda = lambda, 
           lambda.discrete = FALSE, 
           lambda.discrete.num = 100, 
           max.bb.eval = 140,
           min.mesh.size.integer = "1", 
           min.mesh.size.real = paste("r",sqrt(.Machine$double.eps),sep=""),
           min.poll.size.integer = "1", 
           min.poll.size.real = "1",
           nmulti = 0,
           opts=list(),
           random.seed = 42,
           segments = segments, 
           segments.max = 10, 
           segments.min = 1, 
           singular.ok = FALSE,
           tau = NULL,
           weights = NULL)
}

\arguments{

  \item{xz}{ continuous and/or nominal/ordinal
    (\code{\link{factor}}/\code{\link{ordered}}) predictors }

  \item{y}{
    continuous univariate vector
  }

  \item{basis}{ a character string (default \code{basis="additive"})
    indicating whether the additive or tensor product B-spline basis
    matrix for a multivariate polynomial spline or generalized B-spline
    polynomial basis should be used. Note this can be automatically
    determined by cross-validation if \code{cv=TRUE} and
    \code{basis="auto"}, and is an \sQuote{all or none} proposition
    (i.e. interaction terms for all predictors or for no predictors
    given the nature of \sQuote{tensor products}). Note also that if
    there is only one predictor this defaults to \code{basis="additive"}
    to avoid unnecessary computation as the spline bases are equivalent
    in this case }

  \item{complexity}{a character string (default
  \code{complexity="degree-knots"}) indicating whether model
  \sQuote{complexity} is determined by the degree of the spline or by
  the number of segments (\sQuote{knots}). This option allows the user
  to use cross-validation to select either the spline degree (number
  of knots held fixed) or the number of knots (spline degree held
  fixed) or both the spline degree and number of knots}

  \item{cv.df.min}{ the minimum degrees of freedom to allow when
    conducting cross-validation (default \code{cv.df.min=1})}

  \item{cv.func}{a character string (default \code{cv.func="cv.ls"})
    indicating which method to use to select smoothing
    parameters. \code{cv.gcv} specifies generalized cross-validation
    (Craven and Wahba (1979)), \code{cv.aic} specifies expected
    Kullback-Leibler cross-validation (Hurvich, Simonoff, and Tsai
    (1998)), and \code{cv.ls} specifies least-squares
    cross-validation  }
  
  \item{degree}{ integer/vector specifying the degree of the B-spline
  basis for each dimension of the continuous \code{x}}

  \item{degree.max}{ the maximum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.max=10})}
  
  \item{degree.min}{ the minimum degree of the B-spline basis for
    each of the continuous predictors (default \code{degree.min=0})}

  \item{display.nomad.progress}{ a logical value indicating whether to
  display the progress of the NOMAD solver (default \code{display.nomad.progress=TRUE}) }

  \item{display.warnings}{ a logical value indicating whether to
  display warnings (default \code{display.warnings=TRUE}) }

  \item{initial.mesh.size.integer}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{initial.mesh.size.real}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

   \item{knots}{ a character string (default \code{knots="quantiles"})
  specifying where knots are to be placed. \sQuote{quantiles} specifies
  knots placed at equally spaced quantiles (equal number of observations
  lie in each segment) and \sQuote{uniform} specifies knots placed at
  equally spaced intervals. If \code{knots="auto"}, the knot type will
  be automatically determined by cross-validation }

  \item{lambda}{ real/vector for the categorical predictors. If it is
  not NULL, it will be the starting value(s) for lambda}

  \item{lambda.discrete}{ if \code{lambda.discrete=TRUE}, the bandwidth
      will be discretized into \code{lambda.discrete.num+1} points and
      \code{lambda} will be chosen from these points}
  
  \item{lambda.discrete.num}{a positive integer indicating the number of
  discrete values that lambda can assume - this parameter will only be
  used when \code{lambda.discrete=TRUE}}

  \item{max.bb.eval}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details). Default is \code{140} for the aggressive
    speed-oriented profile.
  }

  \item{min.mesh.size.integer}{
    arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{min.mesh.size.real}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{min.poll.size.integer}{
    arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }      

  \item{min.poll.size.real}{
    arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points (default \code{nmulti=0})
  }

  \item{opts}{ list of optional arguments to be passed to
  \code{\link{snomadr}}. If not user-specified, this function applies
  the NOMAD4 path defaults \code{QUAD_MODEL_SEARCH="no"},
  \code{EVAL_QUEUE_SORT="DIR_LAST_SUCCESS"}, and
  \code{DIRECTION_TYPE="ORTHO 2N"} for an aggressive, speed-oriented
  mixed-integer search profile in this specific \code{krscvNOMAD}
  path. User-supplied \code{opts} entries always take precedence. }

\item{random.seed}{
  when it is not missing and not equal to 0, the initial points  will 
  be generated using this seed when \code{nmulti > 0}
}
      
  \item{segments}{ integer/vector specifying the number of segments of
  the B-spline basis for each dimension of the continuous \code{x}
  (i.e. number of knots minus one)}

  \item{segments.max}{ the maximum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.max=10})}
  
  \item{segments.min}{ the minimum segments of the B-spline basis for
    each of the continuous predictors (default \code{segments.min=1})}

  \item{singular.ok}{
    a logical value (default \code{singular.ok=FALSE}) that, when
    \code{FALSE}, discards singular bases during cross-validation (a check
    for ill-conditioned bases is performed).
  }

  \item{tau}{
    if non-null a number in (0,1) denoting the quantile for which a quantile
    regression spline is to be estimated rather than estimating the
    conditional mean (default \code{tau=NULL})
  }
  
  \item{weights}{
    an optional vector of weights to be used in the fitting process.
    Should be \sQuote{NULL} or a numeric vector.  If non-NULL, weighted least
    squares is used with weights \sQuote{weights} (that is, minimizing
    \sQuote{sum(w*e^2)}); otherwise ordinary least squares is used.
  }
}

\details{

  \code{krscvNOMAD} computes NOMAD-based cross-validation for a
  regression spline estimate of a one (1) dimensional dependent variable
  on an \code{r}-dimensional vector of continuous and nominal/ordinal
  (\code{\link{factor}}/\code{\link{ordered}}) predictors.  Numerical
  search for the optimal \code{degree}/\code{segments}/\code{lambda} is
  undertaken using \code{\link{snomadr}}.
  
  The optimal \code{K}/\code{lambda} combination is returned along with
  other results (see below for return values). The method uses kernel
  functions appropriate for categorical (ordinal/nominal) predictors
  which avoids the loss in efficiency associated with sample-splitting
  procedures that are typically used when faced with a mix of continuous
  and nominal/ordinal (\code{\link{factor}}/\code{\link{ordered}})
  predictors.

  For the continuous predictors the regression spline model employs
  either the additive or tensor product B-spline basis matrix for a
  multivariate polynomial spline via the B-spline routines in the GNU
  Scientific Library (\url{https://www.gnu.org/software/gsl/}) and the
  \code{\link{tensor.prod.model.matrix}} function.

  For the discrete predictors the product kernel function is of the
  \sQuote{Li-Racine} type (see Li and Racine (2007) for details).
  
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

  Abramson, M.A. and C. Audet and G. Couture and J.E. Dennis Jr. and
  S. Le Digabel (2011), \dQuote{The NOMAD project}. Software available
  at https://www.gerad.ca/nomad.

  Craven, P. and G. Wahba (1979), \dQuote{Smoothing Noisy Data With
  Spline Functions,} Numerische Mathematik, 13, 377-403.

  Hurvich, C.M. and J.S. Simonoff and C.L. Tsai (1998),
  \dQuote{Smoothing Parameter Selection in Nonparametric Regression
    Using an Improved Akaike Information Criterion,} Journal of the
  Royal Statistical Society B, 60, 271-293.

  Le Digabel, S. (2011), \dQuote{Algorithm 909: NOMAD: Nonlinear
  Optimization With The MADS Algorithm}. ACM Transactions on
  Mathematical Software, 37(4):44:1-44:15.
  
  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
    Theory and Practice,} Princeton University Press.

  Ma, S. and J.S. Racine and L. Yang (2015), \dQuote{Spline
  Regression in the Presence of Categorical Predictors,} Journal of
  Applied Econometrics, Volume 30, 705-717.

  Ma, S. and J.S. Racine (2013), \dQuote{Additive Regression
  Splines with Irrelevant Categorical and Continuous Regressors,}
  Statistica Sinica, Volume 23, 515-541.

}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca} and Zhenghua Nie \email{niez@mcmaster.ca}
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

cv <- krscvNOMAD(x=xdata,y=y,complexity="knots",degree=c(3),segments=c(5))
summary(cv)
}
\keyword{nonparametric}
