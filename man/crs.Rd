\name{crs}
\alias{crs}
\alias{crs.default}
\alias{crs.formula}
\title{Categorical Regression Splines}

\description{ \code{crs} computes a regression spline estimate of a
    one (1) dimensional dependent variable on an \code{r}-dimensional
    vector of continuous and categorical
    (\code{\link{factor}}/\code{\link{ordered}}) predictors (Ma and
    Racine (2013), Ma, Racine and Yang (under revision)).}

\usage{
crs(\dots)
\method{crs}{default}(xz,
    y,
    degree = NULL,
    segments = NULL,
    include = NULL,
    kernel = TRUE,
    lambda = NULL,
    complexity = c("degree-knots","degree","knots"),
    knots = c("quantiles","uniform","auto"),
    basis = c("additive","tensor","glp","auto"),
    deriv = 0,
    data.return = FALSE,
    prune = FALSE,
    model.return = FALSE,
    tau = NULL,
    weights = NULL,
    \dots)

\method{crs}{formula}(formula,
    data = list(),
    degree = NULL,
    segments = NULL,
    include = NULL,
    degree.max = 10, 
    segments.max = 10, 
    degree.min = 0, 
    segments.min = 1, 
    cv.df.min = 1,
    cv = c("nomad","exhaustive","none"),
    cv.threshold = 1000,
    cv.func = c("cv.ls","cv.gcv","cv.aic"),
    kernel = TRUE,
    lambda = NULL,
    lambda.discrete = FALSE,
    lambda.discrete.num = 100,
    complexity = c("degree-knots","degree","knots"),
    knots = c("quantiles","uniform","auto"),
    basis = c("auto","additive","tensor","glp"),
    deriv = 0,
    data.return = FALSE,
    prune = FALSE,
    model.return = FALSE,
    restarts = 0,
    random.seed = 42,
    max.bb.eval = 10000,
    initial.mesh.size.real = "r1.0e-01",
    initial.mesh.size.integer = "1",
    min.mesh.size.real = paste(sqrt(.Machine$double.eps)),
    min.mesh.size.integer = paste(sqrt(.Machine$double.eps)),
    min.poll.size.real = paste(sqrt(.Machine$double.eps)),
    min.poll.size.integer = paste(sqrt(.Machine$double.eps)),
    nmulti = 5,
    tau = NULL,
    weights = NULL,
    singular.ok = FALSE,
    \dots)
}
\arguments{

  \item{xz}{ numeric (\code{x}) and or nominal/ordinal
  (\code{\link{factor}}/\code{\link{ordered}}) predictors (\code{z}) }

\item{y}{ a numeric vector of responses. }

\item{degree}{ integer/vector specifying the polynomial degree of the
B-spline basis for each dimension of the continuous \code{x}
(default \code{degree=3}, i.e. cubic spline)}

\item{segments}{ integer/vector specifying the number of segments of the
B-spline basis for each dimension of the continuous \code{x}
(i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier
curve)}

\item{include}{ integer/vector specifying whether each of the
nominal/ordinal (\code{\link{factor}}/\code{\link{ordered}}) predictors
in \code{x} are included or omitted from the resulting estimate }

  \item{lambda}{ a vector of bandwidths for each dimension of the
    categorical \code{z}}

  \item{lambda.discrete}{ if \code{lambda.discrete=TRUE}, the bandwidth
      will be discretized into \code{lambda.discrete.num+1} points and
      \code{lambda} will be chosen from these points}
  
  \item{lambda.discrete.num}{a positive integer indicating the number of
  discrete values that lambda can assume - this parameter will only be
  used when \code{lambda.discrete=TRUE}}
  
  \item{formula}{ a symbolic description of the model to be fit }
  
  \item{data}{ an optional data frame containing the variables in the
  model }
  
\item{cv}{ a character string (default \code{cv="nomad"}) indicating
  whether to use nonsmooth mesh adaptive direct search, exhaustive
  search, or no search (i.e. use user supplied values for \code{degree},
  \code{segments}, and \code{lambda}) }

\item{cv.threshold}{ an integer (default \code{cv.threshold=1000}) for
  simple problems with no categorical predictors
  (i.e. \code{kernel=FALSE} otherwise
  \code{\link{optim}}/\code{\link{snomadr}} search is necessary) such
  that, if the number of combinations of \code{degree}/\code{segments}
  is less than the threshold and \code{cv="nomad"}, instead use
  exhaustive search (\code{cv="exhaustive"}) }

  \item{cv.func}{a character string (default \code{cv.func="cv.ls"})
    indicating which method to use to select smoothing
    parameters. \code{cv.gcv} specifies generalized cross-validation
    (Craven and Wahba (1979)), \code{cv.aic} specifies expected
    Kullback-Leibler cross-validation (Hurvich, Simonoff, and Tsai
    (1998)), and \code{cv.ls} specifies least-squares
    cross-validation }
  
  \item{kernel}{ a logical value (default \code{kernel=TRUE})
    indicating whether to use kernel smoothing or not }

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

  \item{complexity}{a character string (default
  \code{complexity="degree-knots"}) indicating whether model
  `complexity' is determined by the degree of the spline or by the
  number of segments (i.e. number of knots minus one). This option
  allows the user to use cross-validation to select either the spline
  degree (number of knots held fixed) or the number of knots (spline
  degree held fixed) or both the spline degree and number of knots

  For the continuous predictors the regression spline model employs
  either the additive or tensor product B-spline basis matrix for a
  multivariate polynomial spline via the B-spline routines in the GNU
  Scientific Library (\url{http://www.gnu.org/software/gsl/}) and the
  \code{\link{tensor.prod.model.matrix}} function }

   \item{knots}{ a character string (default \code{knots="quantiles"})
  specifying where knots are to be placed. \sQuote{quantiles} specifies
  knots placed at equally spaced quantiles (equal number of observations
  lie in each segment) and \sQuote{uniform} specifies knots placed at
  equally spaced intervals. If \code{knots="auto"}, the knot type will
  be automatically determined by cross-validation }

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

  \item{deriv}{ an integer \code{l} (default \code{deriv=0}) specifying
  whether to compute the univariate \code{l}th partial derivative for
  each continuous predictor (and difference in levels for each
  categorical predictor) or not and if so what order. Note that if
  \code{deriv} is higher than the spline degree of the associated
  continuous predictor then the derivative will be zero and a warning
  issued to this effect }

  \item{data.return}{ a logical value indicating whether to return
  \code{x,z,y} or not (default \code{data.return=FALSE}) }

  \item{prune}{ a logical value (default \code{prune=FALSE}) specifying
   whether the (final) model is to be \sQuote{pruned} using a stepwise
   cross-validation criterion based upon a modified version of
   \code{\link[MASS]{stepAIC}} (see below for details) }

  \item{model.return}{ a logical value indicating whether to return the
  list of \code{\link{lm}} models or not when \code{kernel=TRUE}
  (default \code{model.return=FALSE}) }

\item{restarts}{ integer specifying the number of times to restart the
    process of finding extrema of the cross-validation function (for the
    bandwidths only) from different (random) initial points }

  \item{random.seed}{ when it is not missing and not equal to 0, the
    initial points will be generated using this seed when using
    \code{\link{frscvNOMAD}} or \code{\link{krscvNOMAD}} and
    \code{nmulti > 0} }
      
  \item{max.bb.eval}{
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

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points (default \code{nmulti=5})
  }
  \item{tau}{
    if non-null a number in (0,1) denoting the quantile for which a quantile
    regression spline is to be estimated rather than estimating the
    conditional mean (default \code{tau=NULL}). Criterion function set
    by \code{cv.func=} are modified accordingly to admit quantile regression.
  }

  \item{weights}{
    an optional vector of weights to be used in the fitting process.
    Should be `NULL' or a numeric vector.  If non-NULL, weighted least
    squares is used with weights `weights' (that is, minimizing
    `sum(w*e^2)'); otherwise ordinary least squares is used.
  }

  \item{singular.ok}{
    a logical value (default \code{singular.ok=FALSE}) that, when
    \code{FALSE}, discards singular bases during cross-validation (a check
    for ill-conditioned bases is performed).
  }
 
  \item{\dots}{ optional arguments }
  
}

\details{

  Typical usages are (see below for a  list of options and also
  the examples at the end of this help file)
  \preformatted{
    ## Estimate the model and let the basis type be determined by
    ## cross-validation (i.e. cross-validation will determine whether to
    ## use the additive, generalized, or tensor product basis)
    
    model <- crs(y~x1+x2)
    
    ## Estimate the model for a specified degree/segment/bandwidth
    ## combination and do not run cross-validation (will use the
    ## additive basis by default)
    
    model <- crs(y~x1+factor(x2),cv="none",degree=3,segments=1,lambda=.1)
    
    ## Plot the mean and (asymptotic) error bounds

    plot(model,mean=TRUE,ci=TRUE)

    ## Plot the first partial derivative and (asymptotic) error bounds

    plot(model,deriv=1,ci=TRUE)    
    }
    
  \code{crs} computes a regression spline estimate of a one (1)
  dimensional dependent variable on an \code{r}-dimensional vector of
  continuous and categorical
  (\code{\link{factor}}/\code{\link{ordered}}) predictors.

  The regression spline model employs the tensor product B-spline basis
  matrix for a multivariate polynomial spline via the B-spline routines
  in the GNU Scientific Library (\url{http://www.gnu.org/software/gsl/})
  and the \code{\link{tensor.prod.model.matrix}} function.

  When \code{basis="additive"} the model becomes additive in nature
  (i.e. no interaction/tensor terms thus semiparametric not fully
  nonparametric).

  When \code{basis="tensor"} the model uses the multivariate tensor
  product basis.

  When \code{kernel=FALSE} the model uses indicator basis functions for
  the nominal/ordinal (\code{\link{factor}}/\code{\link{ordered}})
  predictors rather than kernel weighting.

  When \code{kernel=TRUE} the product kernel function for the discrete
  predictors is of the \sQuote{Li-Racine} type (see Li and Racine (2007)
  for details).

  When \code{cv="nomad"}, numerical search is undertaken using Nonsmooth
  Optimization by Mesh Adaptive Direct Search (Abramson, Audet, Couture,
  Dennis, Jr., and Le Digabel (2011)).
  
  When \code{kernel=TRUE} and \code{cv="exhaustive"}, numerical search
  is undertaken using \code{\link{optim}} and the box-constrained
  \code{L-BFGS-B} method (see \code{\link{optim}} for details). The user
  may restart the algorithm as many times as desired via the
  \code{restarts} argument (default \code{restarts=0}). The approach
  ascends from \code{degree=0} (or \code{segments=0}) through
  \code{degree.max} and for each value of \code{degree} (or
  \code{segments}) searches for the optimal bandwidths. After the most
  complex model has been searched then the optimal
  \code{degree}/\code{segments}/\code{lambda} combination is
  selected. If any element of the optimal \code{degree} (or
  \code{segments}) vector coincides with \code{degree.max} (or \code{segments.max}) a warning
  is produced and the user ought to restart their search with a larger
  value of \code{degree.max} (or \code{segments.max}).

  Note that the default \code{plot} method for a \code{crs} object
  provides some diagnostic measures, in particular, a) residuals versus
  fitted values (useful for checking the assumption that
  \code{E(u|x)=0}), b) a normal quantile-quantile plot which allows
  residuals to be assessed for normality (\code{\link{qqnorm}}), c) a
  scale-location plot that is useful for checking the assumption that
  the errors are iid and, in particular, that the variance is
  homogeneous, and d) \sQuote{Cook's distance} which computes the
  single-case influence function. See below for other arguments for the
  plot function for a \code{crs} object.

  Note that setting \code{prune=TRUE} produces a final \sQuote{pruning}
  of the model via a stepwise cross-validation criterion achieved by
  modifying \code{\link[MASS]{stepAIC}} and replacing \code{extractAIC}
  with \code{extractCV} throughout the function. This option may be
  enabled to remove potentially superfluous bases thereby improving the
  finite-sample efficiency of the resulting model.  Note that if the
  cross-validation score for the pruned model is no better than that for
  the original model then the original model is returned with a warning
  to this effect. Note also that this option can only be used when
  \code{kernel=FALSE}.

}

\value{

  \code{crs} returns a \code{crs} object.  The generic functions
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

  \item{lwr,upr}{ lower/upper bound for a 95\% confidence interval for
  the \code{fitted.values} (conditional mean) obtained from
  \code{\link{predict.lm}} via the argument
  \code{interval="confidence"}}

  \item{residuals}{ residuals computed at the sample points or
  evaluation points }

  \item{degree}{ integer/vector specifying the degree of the B-spline
    basis for each dimension of the continuous \code{x}}
  
  \item{segments}{ integer/vector specifying the number of segments of
  the B-spline basis for each dimension of the continuous \code{x}}

  \item{include}{ integer/vector specifying whether each of the
  nominal/ordinal (\code{\link{factor}}/\code{\link{ordered}})
  predictors \code{z} are included or omitted from the resulting
  estimate if \code{kernel=FALSE} (see below)}

  \item{kernel}{ a logical value indicating whether kernel smoothing was
  used (\code{kernel=TRUE}) or not }
  
  \item{lambda}{ vector of bandwidths used if \code{kernel=TRUE}}
  
  \item{call}{ a symbolic description of the model  }
  
  \item{r.squared}{ coefficient of determination (Doksum and Samarov
  (1995)) }

  \item{model.lm}{ an object of \sQuote{\code{\link{class}}}
  \sQuote{\code{\link{lm}}} if \code{kernel=FALSE} or a list of objects
  of \sQuote{\code{\link{class}}} \sQuote{\code{\link{lm}}} if
  \code{kernel=TRUE} (accessed by \code{model.lm[[1]]},
  \code{model.lm[[2]]},\dots,. By way of example, if \code{foo} is a
  \code{crs} object and \code{kernel=FALSE}, then \code{foo$model.lm}
  is an object of \sQuote{\code{\link{class}}}
  \sQuote{\code{\link{lm}}}, while objects of
  \sQuote{\code{\link{class}}} \sQuote{\code{\link{lm}}} return the
  \code{\link{model.frame}} in \code{model.lm$model} which can be
  accessed via \code{foo$model.lm$model} where \code{foo} is the
  \code{crs} object (the model frame \code{foo$model.lm$model}
  contains the B-spline bases underlying the estimate which might be of
  interest). Again by way of example, when \code{kernel=TRUE} then
  \code{foo$model.lm[[1]]$model} contains the model frame for the first
  unique combination of categorical predictors,
  \code{foo$model.lm[[2]]$model} the second and so forth (the weights
  will potentially differ for each model depending on the value(s) of
  \code{lambda})}

  \item{deriv.mat}{ a matrix of derivatives (or differences in levels
  for the categorical \code{z}) whose order is determined by
  \code{deriv=} in the \code{crs} call }

  \item{deriv.mat.lwr}{ a matrix of 95\% coverage lower bounds for
  \code{deriv.mat}}

  \item{deriv.mat.upr}{ a matrix of 95\% coverage upper bounds for
  \code{deriv.mat}}

  \item{hatvalues}{ the \code{\link{hatvalues}} for the estimated model}

  \item{P.hat}{ the kernel probability estimates corresponding to the
  categorical predictors in the estimated model }

}
  
\section{Usage Issues}{

  Note that when \code{kernel=FALSE} \code{\link{summary}} supports the
  option \code{sigtest=TRUE} that conducts an F-test for significance
  for each predictor.

}
\references{

  Abramson, M.A. and C. Audet and G. Couture and J.E. Dennis Jr. and and
  S. Le Digabel (2011), \dQuote{The NOMAD project}. Software available
  at http://www.gerad.ca/nomad.

  Craven, P. and G. Wahba (1979), \dQuote{Smoothing Noisy Data With
    Spline Functions,} Numerische Mathematik, 13, 377-403.

  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric Estimation of
  Global Functionals and a Measure of the Explanatory Power of
  Covariates in Regression,} The Annals of Statistics, 23 1443-1473.

  Hurvich, C.M. and J.S. Simonoff and C.L. Tsai (1998),
  \dQuote{Smoothing Parameter Selection in Nonparametric Regression
    Using an Improved Akaike Information Criterion,} Journal of the
  Royal Statistical Society B, 60, 271-293.

  Le Digabel, S. (2011), \dQuote{Algorithm 909: NOMAD: Nonlinear
  Optimization With The MADS Algorithm}. ACM Transactions on
  Mathematical Software, 37(4):44:1-44:15.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
    Theory and Practice,} Princeton University Press.

  Ma, S. and J.S. Racine and L. Yang (under revision), \dQuote{Spline
  Regression in the Presence of Categorical Predictors,} Journal of
  Applied Econometrics.

  Ma, S. and J.S. Racine (2013), \dQuote{Additive Regression
  Splines with Irrelevant Categorical and Continuous Regressors,}
  Statistica Sinica, Volume 23, 515-541.

  Racine, J.S. (2011), \dQuote{Cross-Validated Quantile Regression
  Splines,} manuscript.

}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

%\section{Usage Issues}{
%}

\seealso{
  \code{\link{smooth.spline}}, \code{\link{loess}}, \code{\link[np]{npreg}}
}

\examples{
set.seed(42)
## Example - simulated data
n <- 1000
num.eval <- 50
x1 <- runif(n)
x2 <- runif(n)
z <- rbinom(n,1,.5)
dgp <- cos(2*pi*x1)+sin(2*pi*x2)+z
z <- factor(z)
y <- dgp + rnorm(n,sd=.5)

## Estimate a model with specified degree, segments, and bandwidth
model <- crs(y~x1+x2+z,degree=c(5,5),
                       segments=c(1,1),
                       lambda=0.1,
                       cv="none",
                       kernel=TRUE)
summary(model)

## Perspective plot
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(0,num.eval**2),levels=c(0,1)))
z0 <- matrix(predict(model,newdata=newdata),num.eval,num.eval)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2],
                      z=factor(rep(1,num.eval),levels=c(0,1)))
z1 <- matrix(predict(model,newdata=newdata),num.eval,num.eval)
zlim=c(min(z0,z1),max(z0,z1))
persp(x=x1.seq,y=x2.seq,z=z0,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      ticktype="detailed",      
      border="red",
      theta=45,phi=45)
par(new=TRUE)
persp(x=x1.seq,y=x2.seq,z=z1,
      xlab="x1",ylab="x2",zlab="y",zlim=zlim,
      theta=45,phi=45,
      ticktype="detailed",
      border="blue")

## Partial regression surface plot
plot(model,mean=TRUE,ci=TRUE)
\dontrun{
## A plot example where we extract the partial surfaces, confidence
## intervals etc. automatically generated by plot(mean=TRUE,...) but do
## not plot, rather save for separate use.
pdat <- plot(model,mean=TRUE,ci=TRUE,plot.behavior="data")

## Column 1 is the (evaluation) predictor ([,1]), 2-4 ([,-1]) the mean,
## lwr, and upr (note the returned value is a `list' hence pdat[[1]] is
## data for the first predictor, pdat[[2]] the second etc). Note that
## matplot() can plot this nicely.
matplot(pdat[[1]][,1],pdat[[1]][,-1],
        xlab=names(pdat[[1]][1]),ylab=names(pdat[[1]][2]),
        lty=c(1,2,2),col=c(1,2,2),type="l")
}
}
\keyword{nonparametric}

