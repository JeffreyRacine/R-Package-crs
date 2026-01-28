\name{npglpreg}
\alias{npglpreg}
\alias{npglpreg.default}
\alias{npglpreg.formula}

\title{Generalized Local Polynomial Regression}
\description{
  
  \code{npglpreg} computes a generalized local polynomial kernel
    regression estimate (Hall and Racine (2015)) of a one (1)
    dimensional dependent variable on an \code{r}-dimensional vector of
    continuous and categorical
    (\code{\link{factor}}/\code{\link{ordered}}) predictors.
    
  }
\usage{
npglpreg(\dots)

\method{npglpreg}{default}(tydat = NULL,
         txdat = NULL,
         eydat = NULL,
         exdat = NULL,
         bws = NULL,
         degree = NULL,
         leave.one.out = FALSE,
         ckertype = c("gaussian", "epanechnikov", "uniform", "truncated gaussian"),
         ckerorder = 2,
         ukertype = c("liracine", "aitchisonaitken"),
         okertype = c("liracine", "wangvanryzin"),
         bwtype = c("fixed", "generalized_nn", "adaptive_nn", "auto"),
         gradient.vec = NULL,
         gradient.categorical = FALSE,
         cv.shrink = TRUE,
         cv.maxPenalty = sqrt(.Machine$double.xmax),
         cv.warning = FALSE,
         Bernstein = TRUE,
         mpi = FALSE,
         \dots)

\method{npglpreg}{formula}(formula,
         data = list(),
         tydat = NULL,
         txdat = NULL,
         eydat = NULL,
         exdat = NULL,
         bws = NULL,
         degree = NULL,
         leave.one.out = FALSE,
         ckertype = c("gaussian", "epanechnikov","uniform","truncated gaussian"),
         ckerorder = 2,
         ukertype = c("liracine", "aitchisonaitken"),
         okertype = c("liracine", "wangvanryzin"),
         bwtype = c("fixed", "generalized_nn", "adaptive_nn", "auto"),
         cv = c("degree-bandwidth", "bandwidth", "none"),
         cv.func = c("cv.ls", "cv.aic"),
         nmulti = 5,
         random.seed = 42,
         degree.max = 10,
         degree.min = 0,
         bandwidth.max = .Machine$double.xmax,
         bandwidth.min = sqrt(.Machine$double.eps),
         bandwidth.min.numeric = 1.0e-02,
         bandwidth.switch = 1.0e+06,
         bandwidth.scale.categorical = 1.0e+04,
         max.bb.eval = 10000,
         min.epsilon = .Machine$double.eps,
         initial.mesh.size.real = 1,
         initial.mesh.size.integer = 1,
         min.mesh.size.real = sqrt(.Machine$double.eps),
         min.mesh.size.integer = 1, 
         min.poll.size.real = 1,
         min.poll.size.integer = 1, 
         opts=list(),
         restart.from.min = FALSE,
         gradient.vec = NULL,
         gradient.categorical = FALSE,
         cv.shrink = TRUE,
         cv.maxPenalty = sqrt(.Machine$double.xmax),
         cv.warning = FALSE,
         Bernstein = TRUE,
         mpi = FALSE,
         \dots)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit }
  
  \item{data}{ an optional data frame containing the variables in the
    model }

  \item{tydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{txdat}. Defaults to the training data used to
    compute the bandwidth object
  }

  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators. Defaults to the training data used to
    compute the bandwidth object
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    of the dependent variable. Optional, and used only to calculate the
    true errors
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}
  }

  \item{bws}{
    a  vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{txdat}
  }
  
  \item{degree}{ integer/vector specifying the polynomial degree of the
for each dimension of the continuous \code{x} in \code{txdat}}

  \item{leave.one.out}{
    a logical value to specify whether or not to compute the leave one
    out sums. Will not work if \code{exdat} is specified. Defaults to
    \code{FALSE}
  }

  \item{ckertype}{
    character string used to specify the continuous kernel type.
    Can be set as \code{gaussian}, \code{epanechnikov}, or
    \code{uniform}. Defaults to \code{gaussian}.
  }

  \item{ckerorder}{
    numeric value specifying kernel order (one of
    \code{(2,4,6,8)}). Kernel order specified along with a
    \code{uniform} continuous kernel type will be ignored. Defaults to
    \code{2}.
  }

  \item{ukertype}{
    character string used to specify the unordered categorical kernel type.
    Can be set as \code{aitchisonaitken} or \code{liracine}. Defaults to
    \code{liracine}
  }

  \item{okertype}{
    character string used to specify the ordered categorical kernel type.
    Can be set as \code{wangvanryzin} or \code{liracine}. Defaults to
    \code{liracine}
  }

  \item{bwtype}{
    character string used for the continuous variable bandwidth type,
    specifying the type of bandwidth to compute and return in the
    \code{bandwidth} object. If \code{bwtype="auto"}, the bandwidth type
    type will be automatically determined by cross-validation. Defaults
    to \code{fixed}. Option summary:\cr
    \code{fixed}: compute fixed bandwidths \cr
    \code{generalized_nn}: compute generalized nearest neighbor bandwidths \cr
    \code{adaptive_nn}: compute adaptive nearest neighbor bandwidths
  }

  \item{cv}{ a character string (default \code{cv="nomad"}) indicating
    whether to use nonsmooth mesh adaptive direct search, or no search
    (i.e. use supplied values for \code{degree} and \code{bws})
  }

  \item{cv.func}{a character string (default \code{cv.func="cv.ls"})
    indicating which method to use to select smoothing
    parameters. \code{cv.aic} specifies expected Kullback-Leibler
    cross-validation (Hurvich, Simonoff, and Tsai (1998)), and
    \code{cv.ls} specifies least-squares cross-validation
  }

  \item{max.bb.eval}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{min.epsilon}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{initial.mesh.size.real}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{initial.mesh.size.integer}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{min.mesh.size.real}{
    argument passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{min.mesh.size.integer}{
    arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{min.poll.size.real}{
    arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }

  \item{min.poll.size.integer}{
    arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }      

  \item{opts}{ list of optional arguments passed to the NOMAD solver
    (see \code{\link{snomadr}} for further details) }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points (default \code{nmulti=5})
  }
  
  \item{random.seed}{ when it is not missing and not equal to 0, the
    initial points will be generated using this seed when using
    \code{\link{snomadr}} }
      
  \item{degree.max}{ the maximum degree of the polynomial for
    each of the continuous predictors (default \code{degree.max=10})}
  
  \item{degree.min}{ the minimum degree of the polynomial for
    each of the continuous predictors (default \code{degree.min=0})}
  
    \item{bandwidth.max}{ the maximum bandwidth scale (i.e. number of
    scaled standard deviations) for each of the continuous predictors
    (default \code{bandwidth.max=.Machine$double.xmax})}
  
  \item{bandwidth.min}{ the minimum bandwidth scale for each of the
    categorical predictors (default \code{sqrt(.Machine$double.eps)})}

  \item{bandwidth.min.numeric}{ the minimum bandwidth scale (i.e. number
    of scaled standard deviations) for each of the continuous predictors
    (default \code{bandwidth.min=1.0e-02})}

  \item{bandwidth.switch}{ the minimum bandwidth scale (i.e. number of
    scaled standard deviations) for each of the continuous predictors
    (default \code{bandwidth.switch=1.0e+06}) before the local polynomial
    is treated as global during cross-validation at which point a global
    categorical kernel weighted least squares fit is used for
    computational efficiency}

  \item{bandwidth.scale.categorical}{ the upper end for the rescaled
  bandwidths for the categorical predictors (default
  \code{bandwidth.scale.categorical=1.0e+04}) - the aim is to \sQuote{even up}
  the scale of the search parameters as much as possible, so when very
  large scale factors are selected for the continuous predictors, a
  larger value may improve search}

  \item{restart.from.min}{ a logical value indicating to recommence
  \code{\link{snomadr}} with the optimal values found from its first
  invocation (typically quick but sometimes recommended in the field of
  optimization)}
  
  \item{gradient.vec}{ a vector corresponding to the order of the
    partial (or cross-partial) and which variable the partial (or
    cross-partial) derivative(s) are required }

  \item{gradient.categorical}{ a logical value indicating whether
  discrete gradients (i.e. differences in the response from the base
  value for each categorical predictor) are to be computed }

\item{cv.shrink}{ a logical value indicating whether to use ridging
  (Seifert and Gasser (2000)) for ill-conditioned inversion during
  cross-validation (default \code{cv.shrink=TRUE}) or to instead test
  for ill-conditioned matrices and penalize heavily when this is the
  case (much stronger condition imposed on cross-validation)}

\item{cv.maxPenalty}{ a penalty applied during cross-validation when a
  delete-one estimate is not finite or the polynomial basis is not of
  full column rank }

\item{cv.warning}{ a logical value indicating whether to issue an
  immediate warning message when ill conditioned bases are encountered
  during cross-validation (default \code{cv.warning=FALSE}) }

\item{Bernstein}{ a logical value indicating whether to use raw
  polynomials or Bernstein polynomials (default) (note that a Bernstein
  polynomial is also know as a Bezier curve which is also a
  B-spline with no interior knots)}

\item{mpi}{ a logical value (default \code{mpi=FALSE}) that, when
  \code{mpi=TRUE}, can call the \code{npRmpi} rather than the \code{np}
  package (note - code needs to mirror examples in the demo directory of
  the \code{npRmpi} package, you need to broadcast loading of the
  \code{crs} package, and need \code{.Rprofile} in your current
  directory)}

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, training data, and so on, detailed
    below
  }
}
\details{

  Typical usages are (see below for a  list of options and also
  the examples at the end of this help file)
  \preformatted{
    ## Conduct generalized local polynomial estimation
    
    model <- npglpreg(y~x1+x2)
    
    ## Conduct degree 0 local polynomial estimation
    ## (i.e. Nadaraya-Watson)
    
    model <- npglpreg(y~x1+x2,cv="bandwidth",degree=c(0,0))    
    
    ## Conduct degree 1 local polynomial estimation (i.e. local linear)
    
    model <- npglpreg(y~x1+x2,cv="bandwidth",degree=c(1,1))    
    
    ## Conduct degree 2 local polynomial estimation (i.e. local
    ## quadratic)
    
    model <- npglpreg(y~x1+x2,cv="bandwidth",degree=c(2,2))

    ## Plot the mean and bootstrap confidence intervals

    plot(model,ci=TRUE)

    ## Plot the first partial derivatives and bootstrap confidence
    ## intervals

    plot(model,deriv=1,ci=TRUE)

    ## Plot the first second partial derivatives and bootstrap
    ## confidence intervals

    plot(model,deriv=2,ci=TRUE)        
    
  }
  
  This function is in beta status until further notice (eventually it
  will be rolled into the np/npRmpi packages after the final status of
  snomadr/NOMAD gets sorted out).

  Optimizing the cross-validation function jointly for bandwidths
  (vectors of continuous parameters) and polynomial degrees (vectors of
  integer parameters) constitutes a mixed-integer optimization
  problem. These problems are not only \sQuote{hard} from the numerical
  optimization perspective, but are also computationally intensive
  (contrast this to where we conduct, say, local linear regression which
  sets the degree of the polynomial vector to a global value
  \code{degree=1} hence we only need to optimize with respect to the
  continuous bandwidths). Because of this we must be mindful of the
  presence of local optima (the objective function is non-convex and
  non-differentiable). Restarting search from different initial starting
  points is recommended (see \code{nmulti}) and by default this is done
  more than once. We encourage users to adopt \sQuote{multistarting} and
  to investigate the impact of changing default search parameters such
  as \code{initial.mesh.size.real}, \code{initial.mesh.size.integer},
  \code{min.mesh.size.real},
  \code{min.mesh.size.integer},\code{min.poll.size.real}, and
  \code{min.poll.size.integer}. The default values were chosen based on
  extensive simulation experiments and were chosen so as to yield robust
  performance while being mindful of excessive computation - of course,
  no one setting can be globally optimal.

}
\value{
  
  \code{npglpreg} returns a \code{npglpreg} object.  The generic
    functions \code{\link{fitted}} and \code{\link{residuals}} extract
    (or generate) estimated values and residuals. Furthermore, the
    functions \code{\link{summary}}, \code{\link{predict}}, and
    \code{\link{plot}} (options \code{deriv=0}, \code{ci=FALSE}
    [\code{ci=TRUE} produces pointwise bootstrap error bounds],
    \code{persp.rgl=FALSE},
    \code{plot.behavior=c("plot","plot-data","data")},
    \code{plot.errors.boot.num=99},
    \code{plot.errors.type=c("quantiles","standard")}
    [\code{"quantiles"} produces percentiles determined by
    \code{plot.errors.quantiles} below, \code{"standard"} produces error
    bounds given by +/- 1.96 bootstrap standard deviations],
    \code{plot.errors.quantiles=c(.025,.975)}, \code{xtrim=0.0},
    \code{xq=0.5}) support objects of this type. The returned object has
    the following components:

  \item{fitted.values}{ estimates of the regression function
    (conditional mean) at the sample points or evaluation points }

  \item{residuals}{ residuals computed at the sample points or
  evaluation points }

  \item{degree}{ integer/vector specifying the degree of the polynomial
     for each dimension of the continuous \code{x}}
  
  \item{gradient}{ the estimated gradient (vector) corresponding to the vector
    \code{gradient.vec}}

  \item{gradient.categorical.mat}{ the estimated gradient (matrix) for
  the categorical predictors }

  \item{gradient.vec}{ the supplied \code{gradient.vec}}
  
  \item{bws}{ vector of bandwidths }

  \item{bwtype}{ the supplied \code{bwtype}}
  
  \item{call}{ a symbolic description of the model  }
  
  \item{r.squared}{ coefficient of determination (Doksum and Samarov (1995))}

}
\references{

  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric Estimation of
  Global Functionals and a Measure of the Explanatory Power of
  Covariates in Regression,} The Annals of Statistics, 23, 1443-1473.

  Hall, P. and J.S. Racine (2015), \dQuote{Infinite Order
  Cross-Validated Local Polynomial Regression,} Journal of Econometrics,
  185, 510-525.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  Seifert, B. and T. Gasser (2000), \dQuote{Data Adaptive Ridging in
  Local Polynomial Regression,} Journal of Computational and Graphical
  Statistics, 9(2), 338-360.
  

}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca} and Zhenghua Nie \email{niez@mcmaster.ca}
}
\note{
  
  Note that the use of raw polynomials (\code{Bernstein=FALSE}) for
  approximation is appealing as they can be computed and differentiated
  easily, however, they can be unstable (their inversion can be ill
  conditioned) which can cause problems in some instances as the order
  of the polynomial increases. This can hamper search when excessive
  reliance on ridging to overcome ill conditioned inversion becomes
  computationally burdensome.

  \code{npglpreg} tries to detect whether this is an issue or not when
  \code{Bernstein=FALSE} for each \code{numeric} predictor and will
  adjust the search range for \code{\link{snomadr}} and the degree fed
  to \code{npglpreg} if appropriate.

  However, if you suspect that this might be an issue for your specific
  problem and you are using raw polynomials (\code{Bernstein=FALSE}),
  you are encouraged to investigate this by limiting \code{degree.max}
  to value less than the default value (say \code{3}). Alternatively,
  you might consider re-scaling your \code{numeric} predictors to lie in
  \eqn{[0,1]} using \code{\link{scale}}.

  For a given predictor \eqn{x} you can readily determine if this is an
  issue by considering the following: Suppose \eqn{x} is given by

  \preformatted{
    x <- runif(100,10000,11000)
    y <- x + rnorm(100,sd=1000)
  }

  so that a polynomial of order, say, \eqn{5} would be ill
  conditioned. This would be apparent if you considered

  \preformatted{
    X <- poly(x,degree=5,raw=TRUE)
    solve(t(X)\%*\%X)
  }
  
  which will throw an error when the polynomial is ill conditioned,
  or
  
  \preformatted{
    X <- poly(x,degree=5,raw=TRUE)
    lm(y~X)
  }

  which will return \code{NA} for one or more coefficients when this is
  an issue.

  In such cases you might consider transforming your \code{numeric}
  predictors along the lines of the following:
  \preformatted{
    x <- as.numeric(scale(x))
    X <- poly(x,degree=5,raw=TRUE)
    solve(t(X)\%*\%X)
    lm(y~X)    
  }

  Note that now your least squares coefficients (i.e. first derivative
  of \eqn{y} with respect to \eqn{x}) represent the effect of a one
  standard deviation change in \eqn{x} and not a one unit change.

  Alternatively, you can use Bernstein polynomials by not setting
  \code{Bernstein=FALSE}.
  
  
}

\seealso{
\code{\link[np]{npreg}}
}
\examples{
\dontrun{
set.seed(42)
n <- 100
x1 <- runif(n,-2,2)
x2 <- runif(n,-2,2)
y <- x1^3 + rnorm(n,sd=1)

## Ideally the method should choose large bandwidths for x1 and x2 and a
## generalized polynomial that is a cubic for x1 and degree 0 for x2.

model <- npglpreg(y~x1+x2,nmulti=1)
summary(model)

## Plot the partial means and percentile confidence intervals
plot(model,ci=TRUE)
## Extract the data from the plot object and plot it separately
# myplot.dat <- plot(model,plot.behavior="data",ci=TRUE)
# matplot(myplot.dat[[1]][,1],myplot.dat[[1]][,-1],type="l")
# matplot(myplot.dat[[2]][,1],myplot.dat[[2]][,-1],type="l")
} % end dontrun
}
\keyword{ nonparametric }
\keyword{ regression }

