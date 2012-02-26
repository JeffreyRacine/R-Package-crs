\name{crsivderiv}
\alias{crsivderiv}

\title{
Nonparametric Instrumental Regression
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
crsivderiv(y,
      z,
      w,
      x = NULL,
      zeval = NULL,
      weval = NULL,
      xeval = NULL,
      iterate.max = 100,
      iterate.tol = 1.0e-04,
      constant = 0.5,
      start.phi.zero = FALSE,
      stop.on.increase = TRUE,
      opts = list("MAX_BB_EVAL"=10000,
                  "EPSILON"=.Machine$double.eps,
                  "INITIAL_MESH_SIZE"="r1.0e-01",
                  "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                  "MIN_POLL_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
                  "DISPLAY_DEGREE"=0),
      deriv = 0,
      ...)
}

\arguments{
  \item{y}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{z}
  }

  \item{z}{
    a \eqn{p}-variate data frame of endogenous predictors. The data
    types may be continuous, discrete (unordered and ordered factors),
    or some combination thereof
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
  \item{zeval}{
    a \eqn{p}-variate data frame of endogenous predictors on which the
    regression will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{z}
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
\item{iterate.max}{
  an integer indicating the maximum number of iterations permitted
  before termination occurs when using Landweber-Fridman iteration
}
\item{iterate.tol}{
  the search tolerance for the stopping rule when using Landweber-Fridman iteration
}
\item{constant}{
  the constant to use when using  Landweber-Fridman iteration
}
\item{start.phi.zero}{
  a logical value indicating whether to commence Landweber-Fridman
  assuming \eqn{\varphi_{-1}=0}{phi[-1]=0} (proper Landweber-Fridman) or
  instead begin from \eqn{E(y|z)} (defaults to \code{FALSE}, see details
  below)
}
\item{smooth.while.iterating}{
  a logical value indicating whether to re-smooth
  \eqn{E(\varphi|w)}{E(phi|w)} during Landweber-Fridman iteration
  (defaults to \code{TRUE})
}
\item{stop.on.increase}{
  a logical value (defaults to \code{TRUE}) indicating whether to halt
  iteration if the stopping criterion (see below) increases over the
  course of one iteration (i.e. it may be above the iteration tolerance
  but increased)
}
\item{opts}{
  arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
  further details)
}
\item{deriv}{ an integer \code{l} (default \code{deriv=0}) specifying
  whether to compute the univariate \code{l}th partial derivative for
  each continuous predictor (and difference in levels for each
  categorical predictor) or not and if so what order. Note that if
  \code{deriv} is higher than the spline degree of the associated
  continuous predictor then the derivative will be zero and a warning
  issued to this effect
}
\item{...}{
  additional arguments supplied to \code{\link{crs}}
}
}
\details{

  For Landweber-Fridman iteration, an optimal stopping rule based upon
  \eqn{||E(y|w)-E(\varphi_k(z,x)|w)||^2 }{||E(y|w)-E(phi(z,x)|w)||^2 }
  is used to terminate iteration. However, if local rather than global
  optima are encountered the resulting estimates can be overly noisy. To
  best guard against this eventuality set \code{nmulti} to a larger
  number than the default \code{nmulti=5} for \code{\link{crs}} when
  using \code{cv="nomad"} or instead use \code{cv="exhaustive"} if
  possible (this may not be feasible for non-trivial problems).

  When using Landweber-Fridman iteration, iteration will terminate
  when either the change in the value of
  \eqn{||(E(y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } from iteration to iteration is
  less than \code{iterate.tol} or we hit \code{iterate.max} or
  \eqn{||E((y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } stops falling in value and
  starts rising.

  When your problem is a simple one (e.g. univariate \eqn{Z}, \eqn{W},
  and \eqn{X}) you might want to avoid \code{cv="nomad"} and instead use
  \code{cv="exhaustive"} since exhaustive search may be feasible (for
  \code{degree.max} and \code{segments.max} not overly large). This will
  guarantee an exact solution for each iteration (i.e. there will be no
  errors arising due to numerical search).

  \code{demo(crsivderiv)}, \code{demo(crsivderiv_exog)}, and
  \code{demo(crsivderiv_exog_persp)} provide flexible interactive
  demonstrations similar to the example below that allow you to modify
  and experiment with parameters such as the sample size, method, and so
  forth in an interactive session.
  
}
\value{
  
  \code{crsivderiv} returns a \code{\link{crs}} object.  The generic
    functions \code{\link{fitted}} and \code{\link{residuals}} extract
    (or generate) estimated values and residuals. Furthermore, the
    functions \code{\link{summary}}, \code{\link{predict}}, and
    \code{\link{plot}} (options \code{mean=FALSE}, \code{deriv=i} where
    \eqn{i} is an integer, \code{ci=FALSE},
    \code{plot.behavior=c("plot","plot-data","data")}) support objects
    of this type.
    
    See \code{\link{crs}} for details on the return object components.

    In addition to the standard \code{\link{crs}} components,
    \code{crsivderiv} returns components \code{phi.prime}, \code{phi},
    \code{num.iterations} and \code{norm.stop}.
    
}
  
\references{
  
  Carrasco, M. and J.P. Florens and E. Renault (2007), \dQuote{Linear
  Inverse Problems in Structural Econometrics Estimation Based on
  Spectral Decomposition and Regularization,} In: James J. Heckman and
  Edward E. Leamer, Editor(s), Handbook of Econometrics, Elsevier, 2007,
  Volume 6, Part 2, Chapter 77, Pages 5633-5751

  Darolles, S. and Y. Fan and J.P. Florens and E. Renault (forthcoming),
  \dQuote{Nonparametric Instrumental Regression,} Econometrica.

  Feve, F. and J.P. Florens (2010), \dQuote{The practice of
  non-parametric estimation by solving inverse problems: the example of
  transformation models,} Econometrics Journal, 13, S1-S27.

  Florens, J.P. and J.S. Racine (2012), \dQuote{Nonparametric
  Instrumental Variables Estimation of Marginal Effects,} Working Paper.

  Fridman, V. M. (1956), \dQuote{A Method of Successive Approximations
  for Fredholm Integral Equations of the First Kind,} Uspeskhi,
  Math. Nauk., 11, 233-334, in Russian.

  Horowitz, J.L. (2011), \dQuote{Applied Nonparametric Instrumental
  Variables Estimation,}, Econometrica, 79, 347-394.
  
  Landweber, L. (1951), \dQuote{An iterative formula for Fredholm
  integral equations of the first kind,} American Journal of
  Mathematics, 73, 615-24.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}, Samuele Centorrino
 \email{samuele.centorrino@univ-tlse1.fr}
}
\note{
This function should be considered to be in `beta test' status until further notice.
}

\seealso{
\code{\link[np]{npreg}}, \code{\link{crsiv}}, \code{\link{crs}}
}
\examples{
## This illustration was made possible by Samuele Centorrino
## <samuele.centorrino@univ-tlse1.fr>

set.seed(42)
n <- 1500

## Turn off display I/O for crs() NOMAD search (i.e. keep the standard
## defaults and add "DISPLAY_DEGREE"=0)...

opts <- list("MAX_BB_EVAL"=10000,
             "EPSILON"=.Machine$double.eps,
             "INITIAL_MESH_SIZE"="r1.0e-01",
             "MIN_MESH_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
             "MIN_POLL_SIZE"=paste("r",sqrt(.Machine$double.eps),sep=""),
             "DISPLAY_DEGREE"=0)

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

model.ivderiv <- crsivderiv(y=y,z=z,w=w)

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
}
\keyword{ instrument }
