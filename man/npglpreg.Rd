\name{npglpreg}
\alias{npglpreg}
\alias{npglpreg.default}
\alias{npglpreg.formula}

\title{Generalized Local Polynomial Regression}
\description{
  
  \code{npglpreg} computes a generalized local polynomial kernel
    regression estimate (Hall and Racine (2011)) of a one (1)
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
         ukertype = c("liracine", "aitchisonaitken"),
         okertype = c("liracine", "wangvanryzin"),
         bwtype = c("fixed", "generalized_nn", "adaptive_nn"),
         raw = TRUE,
         gradient.vec = NULL,
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
         ukertype = c("liracine", "aitchisonaitken"),
         okertype = c("liracine", "wangvanryzin"),
         bwtype = c("fixed", "generalized_nn", "adaptive_nn"),
         cv = c("degree-bandwidth", "bandwidth", "none"),
         cv.func = c("cv.ls", "cv.gcv", "cv.aic"),
         opts = list(MAX_BB_EVAL = 10000,
           EPSILON = .Machine$double.eps,
           INITIAL_MESH_SIZE = "1.0e-01",
           MIN_MESH_SIZE = paste("r", sqrt(.Machine$double.eps), sep = ""),
           MIN_POLL_SIZE = paste("r", sqrt(.Machine$double.eps), sep = "")),
         nmulti = 5,
         random.seed = 42,
         degree.max = 5,
         degree.min = 0,
         bandwidth.max = 1.0e+05,
         bandwidth.min = 1.0e-03,
         raw = TRUE,
         gradient.vec = NULL,
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
    compute the bandwidth object.
  }

  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    of the dependent variable. Optional, and used only to calculate the
    true errors. 
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{bws}{
    a  vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{txdat}.
  }
  
  \item{degree}{ integer/vector specifying the polynomial degree of the
for each dimension of the continuous \code{x} in \code{txdat}}

  \item{leave.one.out}{
    a logical value to specify whether or not to compute the leave one
    out sums. Will not work if \code{exdat} is specified. Defaults to
    \code{FALSE}.
  }

  \item{ukertype}{
    character string used to specify the unordered categorical kernel type.
    Can be set as \code{aitchisonaitken} or \code{liracine}. Defaults to
    \code{liracine}.
  }

  \item{okertype}{
    character string used to specify the ordered categorical kernel type.
    Can be set as \code{wangvanryzin} or \code{liracine}. Defaults to
    \code{liracine}.
  }

  \item{bwtype}{
    character string used for the continuous variable bandwidth type,
    specifying the type of bandwidth to compute and return in the
    \code{bandwidth} object. Defaults to \code{fixed}. Option
    summary:\cr
    \code{fixed}: compute fixed bandwidths \cr
    \code{generalized_nn}: compute generalized nearest neighbors \cr
    \code{adaptive_nn}: compute adaptive nearest neighbors
    
  }
  \item{cv}{
  a character string (default \code{cv="nomad"}) indicating whether to
  use nonsmooth mesh adaptive direct search, or no search (i.e. use 
  supplied values for \code{degree} and \code{bws})
}
  \item{cv.func}{a character string (default \code{cv.func="cv.ls"})
    indicating which method to use to select smoothing
    parameters. \code{cv.gcv} specifies generalized cross-validation
    (Craven and Wahba (1979)), \code{cv.aic} specifies expected
    Kullback-Leibler cross-validation (Hurvich, Simonoff, and Tsai
    (1998)), and \code{cv.ls} specifies least-squares
    cross-validation }
  \item{opts}{
    arguments passed to the NOMAD solver (see \code{\link{snomadr}} for
    further details)
  }
  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points (default \code{nmulti=5})
  }
  \item{random.seed}{ when it is not missing and not equal to 0, the
    initial points will be generated using this seed when using
    \code{\link{snomadr}} }
      
  \item{degree.max}{ the maximum degree of the polynomial for
    each of the continuous predictors (default \code{degree.max=5})}
  
  \item{degree.min}{ the minimum degree of the polynomial for
    each of the continuous predictors (default \code{degree.min=0})}
  
    \item{bandwidth.max}{ the maximum bandwidth scale (i.e. number of
    standard deviations) for each of the continuous predictors (default
    \code{bandwidth.max=10000})}
  
  \item{bandwidth.min}{ the minimum bandwidth scale (i.e. number of
    standard deviations) for each of the continuous predictors (default
    \code{bandwidth.min=0.001})}
  
  \item{raw}{
  when TRUE evaluate raw polynomials of degree 1 to \code{degree} over the
  specified set of points \code{x} otherwise evaluate orthogonal polynomials
  (orthogonal to the constant polynomial of degree 0)
}
\item{gradient.vec}{
  a vector corresponding to the order of the partial (or cross-partial)
  and which variable the partial (or cross-partial) derivative(s) are required
}
  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, training data, and so on, detailed
    below. 
  }
}
\details{
  This function is in beta status until further notice.

  Important note - prediction with newdata only works if there exists
  more than one row in the newdata data frame.
}
\value{
  \code{npglpreg} returns a \code{npglpreg} object.  The generic
    functions \code{\link{fitted}} and \code{\link{residuals}} extract
    (or generate) estimated values and residuals. Furthermore, the
    functions \code{\link{summary}}, \code{\link{predict}}, and
    \code{\link{plot}} (options \code{deriv=0}, \code{ci=FALSE},
    \code{persp.rgl=FALSE},
    \code{plot.behavior=c("plot","plot-data","data")}) support objects
    of this type. The returned object has the following components:

  \item{fitted.values}{ estimates of the regression function
    (conditional mean) at the sample points or evaluation points }

  \item{residuals}{ residuals computed at the sample points or
  evaluation points }

  \item{degree}{ integer/vector specifying the degree of the polynomial
     for each dimension of the continuous \code{x}}
  
  \item{kernel}{ a logical value indicating whether kernel smoothing was
    used (\code{kernel=TRUE}) or not }

  \item{gradient}{ the estimated gradient corresponding to the vector \code{gradient.vec}}

  \item{gradient.vec}{ the supplied \code{gradient.vec}}
  
  \item{bws}{ vector of bandwidths }

  \item{bwtype}{ the supplied \code{bwtype}}
  
  \item{call}{ a symbolic description of the model  }
  
  \item{r.squared}{ coefficient of determination }

}
\references{
  
  Hall, P. and J.S. Racine (2011), \dQuote{Cross-Validated
  Generalized Local Polynomial Regression,} Manuscript.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca} and Zhenghua Nie <niez@mcmaster.ca>
}
\note{
TBA
}

\seealso{
\code{\link[np]{npreg}}
}
\examples{
set.seed(42)
n <- 100
x1 <- runif(n)
x2 <- runif(n)
y <- x1^3 + rnorm(n,sd=.1)

## Ideally the method should choose large bandwidths for x1 and x2 and a
## generalized polynomial that is a cubic for x1 and degree 0 for x2.

model <- npglpreg(y~x1+x2)
summary(model)
}
\keyword{ nonparametric }
\keyword{ regression }

