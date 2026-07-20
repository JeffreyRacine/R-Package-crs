\name{crsivderiv}
\alias{crsivderiv}
\alias{crsivderiv.default}

\title{
Nonparametric Instrumental Derivatives
}
\description{
  
  \code{crsivderiv} uses the approach of Florens and Racine (2012) to
  compute the partial derivative of a nonparametric estimation of an
  instrumental regression function \eqn{\varphi}{phi} defined by
  conditional moment restrictions stemming from a structural econometric
  model: \eqn{E [Y - \varphi (Z,X) | W ] = 0}{E [Y - phi (Z,X) | W ] =
  0}, and involving endogenous variables \eqn{Y} and \eqn{Z} and
  exogenous variables \eqn{X} and instruments \eqn{W}. The derivative
  function \eqn{\varphi'}{phi'} is the solution of an ill-posed inverse
  problem, and is computed using Landweber-Fridman regularization.

}
\usage{
crsivderiv(y, \dots)

\method{crsivderiv}{default}(y,
           z,
           w,
           x = NULL,
           zeval = NULL,
           weval = NULL,
           xeval = NULL,
           constant = 0.5,
           display.nomad.progress = TRUE,
           display.warnings = TRUE,
           iterate.diff.tol = 1.0e-08,
           iterate.max = 1000,
           opts = list("MAX_BB_EVAL"=10000,
                       "EPSILON"=.Machine$double.eps,
                       "INITIAL_MESH_SIZE"="r1.0e-01",
                       "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                       "MIN_FRAME_SIZE"=paste("r",1,sep=""),
                       "DISPLAY_DEGREE"=0),
           penalize.iteration = TRUE,
           smooth.residuals = TRUE,
           start.from = c("Eyz","EEywz"),
           starting.values = NULL,
           stop.on.increase = TRUE,
           \dots)
}

\arguments{
\subsection{Data, Model Inputs And Formula Interface}{
These arguments identify the response, endogenous variables, instruments, and exogenous covariates.
}

  \item{w}{
    a \eqn{q}-variate data frame of instruments. The data types may be
    continuous, discrete (unordered and ordered factors), or some
    combination thereof
  }

  \item{x}{
    an \eqn{r}-variate data frame of exogenous predictors. The data
    types may be continuous, discrete (unordered and ordered factors),
    or some combination thereof
  }


  \item{y}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{z}
  }


  \item{z}{
    a one-column data frame of continuous endogenous predictors. The
    current implementation of \code{crsivderiv} supports univariate
    continuous \eqn{z} only
  }

\subsection{Evaluation Inputs}{
These arguments identify evaluation data for the derivative fit.
}

  \item{weval}{
    a \eqn{q}-variate data frame of instruments on which the regression
    will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{w}
  }

  \item{xeval}{
    an \eqn{r}-variate data frame of exogenous predictors on which the
    regression will be estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{x}
  }


  \item{zeval}{
    a one-column data frame of continuous endogenous predictors on
    which the regression will be estimated (evaluation data). By
    default, evaluation takes place on the data provided by \code{z}
  }

\subsection{Landweber-Fridman Iteration Controls}{
These arguments control iteration, residual smoothing, starting values, and stopping behavior.
}

\item{constant}{
  the constant to use when using  Landweber-Fridman iteration
}

\item{iterate.diff.tol}{
  the search tolerance for the difference in the stopping rule from
  iteration to iteration when using Landweber-Fridman
  (disable by setting to zero)
}

\item{iterate.max}{
  an integer indicating the maximum number of complete positive-index
  Landweber-Fridman states that may be evaluated. The initialization state
  \eqn{N=0} is not counted
}

\item{penalize.iteration}{ a logical value indicating whether to
    penalize the norm by the number of iterations or not (default
    \code{TRUE})
}    

\item{smooth.residuals}{
  a logical value (defaults to \code{TRUE}) indicating whether to
  optimize bandwidths for the regression of \eqn{y-\varphi(z)}{y-phi(z)}
  on \eqn{w} or for the regression of \eqn{\varphi(z)}{phi(z)} on
  \eqn{w} during iteration
 }

\item{start.from}{ a character string indicating whether to start from
  \eqn{E(Y|z)} (default, \code{"Eyz"}) or from \eqn{E(E(Y|z)|z)} (this can
  be overridden by providing \code{starting.values} below)
}  

\item{starting.values}{ optional derivative values for the initialization
  state \eqn{\varphi'_0}{phi'[0]}. They are used directly before the first
  reported update; no hidden \eqn{-1} to \eqn{0} update is performed. When
  \code{NULL}, the initialization is obtained from \eqn{E(y|z)} according to
  \code{start.from} (see details below)
}

\item{stop.on.increase}{
  a logical value (defaults to \code{TRUE}) indicating whether to halt
  iteration if the stopping criterion (see below) increases over the
  course of one iteration (i.e. it may be above the iteration tolerance
  but increased)
}


\subsection{Warnings And Progress}{
These arguments control warnings and displayed optimizer progress.
}

\item{display.nomad.progress}{ a logical value indicating whether to
  display the progress of the NOMAD solver (default \code{display.nomad.progress=TRUE}) }


\item{display.warnings}{ a logical value indicating whether to
  display warnings (default \code{display.warnings=TRUE}) }
\subsection{Additional Arguments}{
Further NOMAD and CRS controls are passed through to lower-level routines.
}

\item{...}{
  additional arguments supplied to \code{\link{crs}}
}


\item{opts}{
  arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
  further details)
}

}
\details{

  The iteration index has a single meaning throughout the returned object.
  The initialization consists of \eqn{\varphi'_0}{phi'[0]} and its centered
  integral \eqn{\varphi_0}{phi[0]}. For \eqn{N=1,2,\ldots}, one adjoint
  update produces \eqn{\varphi'_N}{phi'[N]}, integration and centering produce
  \eqn{\varphi_N}{phi[N]}, and the residual at that same curve produces the
  stopping value for state \eqn{N}. Thus column \eqn{N} of
  \code{phi.prime.mat} and \code{phi.mat}, and element \eqn{N} of
  \code{norm.stop}, always describe one coherent state. With
  \code{penalize.iteration=TRUE}, that element is \eqn{N q_N}, where
  \eqn{q_N} is the normalized squared residual criterion at
  \eqn{\varphi_N}{phi[N]}.

  For Landweber-Fridman iteration, an optimal stopping rule based upon
  \eqn{||E(y|w)-E(\varphi_k(z,x)|w)||^2 }{||E(y|w)-E(phi(z,x)|w)||^2 }
  is used to terminate iteration. However, if local rather than global
  optima are encountered the resulting estimates can be overly noisy. To
  best guard against this eventuality set \code{nmulti} to a larger
  number than the default \code{nmulti=2} for \code{\link{crs}} when
  using \code{cv="nomad"} or instead use \code{cv="exhaustive"} if
  possible (this may not be feasible for non-trivial problems).

  Note that for subsequent Landweber-Fridman iterations, a \dQuote{warm
  start} strategy is employed. The optimal parameters (spline degree,
  number of segments, and bandwidths or inclusion indicators) from the
  previous iteration are used as starting values for the current
  iteration. The user-supplied \code{nmulti} is respected for all
  iterations. For iterations after the first successful one, these
  optimal parameters serve as the first of the multiple initial points
  (a warm start), while any remaining restarts are cold starts. If
  \code{nmulti} is not explicitly supplied by the user, it defaults to
  the \code{\link{crs}} default (2) for the first iteration and to 1 for
  all subsequent iterations. This strategy provides a balance between
  computational efficiency and robustness, allowing the NOMAD solver to
  refine the structural parameters as the residuals evolve
  incrementally while still guarding against local optima.

  When using Landweber-Fridman iteration, iteration will terminate
  when either the change in the value of
  \eqn{||(E(y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } from iteration to iteration is
  less than \code{iterate.diff.tol} or we hit \code{iterate.max} or
  \eqn{||(E(y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } stops falling in value and
  starts rising.

  When your problem is a simple one (e.g. univariate \eqn{Z}, \eqn{W},
  and \eqn{X}) you might want to avoid \code{cv="nomad"} and instead use
  \code{cv="exhaustive"} since exhaustive search may be feasible (for
  \code{degree.max} and \code{segments.max} not overly large). This will
  guarantee an exact solution for each iteration (i.e. there will be no
  errors arising due to numerical search).

  The current implementation supports a single continuous endogenous
  regressor only. Instrument and exogenous regressor data may still be
  mixed continuous and categorical.

}
\value{
  
  \code{crsivderiv} returns a \code{crsivderiv} object (which inherits
    from the \code{\link{crs}} class).  The generic functions
    \code{\link{print}}, \code{\link{summary}}, \code{\link{fitted}},
    \code{\link{residuals}}, \code{\link{predict}}, and
    \code{\link{plot}} support objects of this type.

    For the \code{\link{plot}} function, the options are
    \code{plot.data=FALSE} or \code{data_overlay=FALSE} (logical values
    indicating whether to plot the data as a scatter plot),
    \code{phi=FALSE} (a logical value indicating whether to plot the
    reconstructed structural function rather than its derivative), and
    \code{output=c("plot","data","plot-data")} or \code{behavior}
    (whether to draw, return the plot data, or both). Bootstrap,
    asymptotic interval, surface-rendering, rug, legend, and regression
    gradient plot controls are not supported for this curve route and
    fail fast when supplied. See \code{\link{plot.crs}} for the shared CRS
    plot-output and data-overlay conventions. Note that the \code{\link{plot}}
    method for \code{crsivderiv} objects currently only supports a univariate
    continuous endogenous predictor \eqn{z}.
    
    See \code{\link{crs}} for details on the return object components.

    In addition to the standard \code{\link{crs}} components,
    \code{crsivderiv} returns components \code{phi.prime}, \code{phi},
    \code{phi.prime.mat}, \code{phi.mat}, \code{num.iterations},
    \code{norm.stop}, \code{norm.value} and \code{convergence}.
    \code{num.iterations} is the selected positive-index state, and
    \code{phi.prime} and \code{phi} are taken from that same column of their
    respective matrices. The matrices and \code{norm.stop} contain all
    evaluated states, so their common length can exceed
    \code{num.iterations} when a later evaluated state triggers stopping or
    the stopping-rule selector chooses an earlier state.
    
}
  
\references{
  
  Carrasco, M. and J.P. Florens and E. Renault (2007), \dQuote{Linear
  Inverse Problems in Structural Econometrics Estimation Based on
  Spectral Decomposition and Regularization,} In: James J. Heckman and
  Edward E. Leamer, Editor(s), Handbook of Econometrics, Elsevier, 2007,
  Volume 6, Part 2, Chapter 77, Pages 5633-5751

  Darolles, S. and Y. Fan and J.P. Florens and E. Renault (2011),
  \dQuote{Nonparametric Instrumental Regression,} Econometrica, 79,
  1541-1565.
  
  Feve, F. and J.P. Florens (2010), \dQuote{The Practice of
  Non-parametric Estimation by Solving Inverse Problems: The Example of
  Transformation Models,} Econometrics Journal, 13, S1-S27.

  Florens, J.P. and J.S. Racine (2012), \dQuote{Nonparametric
  Instrumental Derivatives,} Working Paper.

  Fridman, V. M. (1956), \dQuote{A Method of Successive Approximations
  for Fredholm Integral Equations of the First Kind,} Uspeskhi,
  Math. Nauk., 11, 233-334, in Russian.

  Horowitz, J.L. (2011), \dQuote{Applied Nonparametric Instrumental
  Variables Estimation,} Econometrica, 79, 347-394.
  
  Landweber, L. (1951), \dQuote{An Iterative Formula for Fredholm
  Integral Equations of the First Kind,} American Journal of
  Mathematics, 73, 615-24.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\note{
This function currently supports univariate \code{z} only.
This function should be considered to be in \sQuote{beta test} status until
further notice.
}

\seealso{
\code{npreg}, \code{\link{crsiv}}, \code{\link{crs}},
\code{\link{plot.crs}}
}
\examples{
\dontrun{
## This illustration was made possible by Samuele Centorrino
## <samuele.centorrino@univ-tlse1.fr>

set.seed(42)
n <- 1000

## For trimming the plot (trim .5\% from each tail)

trim <- 0.005

## The DGP is as follows:

## 1) y = phi(z) + u

## 2) E(u|z) != 0 (endogeneity present)

## 3) Suppose there exists an instrument w such that z = f(w) + v and
## E(u|w) = 0

## 4) We generate v, w, and generate u such that u and z are
## correlated. To achieve this we express u as a function of v (i.e. u =
## gamma v + eps)

v <- rnorm(n,mean=0,sd=0.27)
eps <- rnorm(n,mean=0,sd=0.05)
u <- -0.5*v + eps
w <- rnorm(n,mean=0,sd=1)

## In Darolles et al (2011) there exist two DGPs. The first is
## phi(z)=z^2 and the second is phi(z)=exp(-abs(z)) (which is
## discontinuous and has a kink at zero).

fun1 <- function(z) { z^2 }
fun2 <- function(z) { exp(-abs(z)) }

z <- 0.2*w + v

## Generate two y vectors for each function.

y1 <- fun1(z) + u
y2 <- fun2(z) + u

## You set y to be either y1 or y2 (ditto for phi) depending on which
## DGP you are considering:

y <- y1
phi <- fun1

## Sort on z (for plotting)

ivdata <- data.frame(y,z,w,u,v)
ivdata <- ivdata[order(ivdata$z),]
rm(y,z,w,u,v)
attach(ivdata)

## Setting cv.threshold = 0 forces NOMAD search instead of exhaustive search
## when no categorical predictors are present. This avoids unnecessary
## evaluation of all degree/segment combinations in the examples and, for
## crsiv() and crsivderiv(), ensures that the warm-start strategy is used.
model.ivderiv <- crsivderiv(y=y,z=z,w=w,cv.threshold=0)

ylim <-c(quantile(model.ivderiv$phi.prime,trim),
         quantile(model.ivderiv$phi.prime,1-trim))

plot(z,model.ivderiv$phi.prime,
     xlim=quantile(z,c(trim,1-trim)),
     main="",
     ylim=ylim,
     xlab="Z",
     ylab="Derivative",
     type="l",
     lwd=2)
rug(z)
} % end dontrun
}
\keyword{ instrument }
