% $ID$
\name{frscv}
\alias{frscv}

\title{Categorical Factor Regression Spline Cross-Validation}

\description{
  
  \code{frscv} computes cross-validation for a regression spline
  estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous predictors and
  ordinal/nominal factor predictors.
  
}

\usage{
frscv(xz,
      y,
      basis.maxdim = 10,
      complexity = c("degree","knots"),
      basis = c("additive-tensor","additive","tensor","auto"),
      cv.norm = c("L2","L1"),
      degree = 3,
      nbreak = 1)
}

\arguments{

  \item{y}{
    continuous univariate vector
  }
  
  \item{xz}{
    continuous and/or ordinal/nominal predictors
  }

  \item{basis.maxdim}{
    integer specifying the maximum degree of the regression spline,
    defaults to \code{10}
  }

  \item{complexity}{a character string (default \code{complexity="degree"})
  indicating whether model `complexity' is determined by the degree of
  the spline or by the number of breakpoints (`knots'). This option
  allows the user to use cross-validation to select either the spline
  degree (number of knots held fixed) or the number of knots (spline
  degree held fixed)}

  \item{basis}{ a character string (default \code{basis="additive-tensor"})
  indicating whether tensor product terms (i.e. interaction terms) are
  included in the model or not. Note this is automatically determined
  when \code{basis="auto"}, and is an \sQuote{all or none}
  proposition (i.e. interaction terms for all predictors or for no
  predictors) }

  \item{cv.norm}{ a character string (default \code{cv.norm="L2"}
  indicating whether the L2-norm is used (sum of squared delete-one
  residuals) or the L1-norm (sum of absolute delete-one-residuals) in
  the cross-validation function }

  \item{degree}{ integer/vector specifying the degree of the regression
    spline for each dimension of the continuous \code{x}}
  
  \item{nbreak}{ integer/vector specifying the number of breakpoints of
  the regression spline for each dimension of the continuous \code{x}}
  
}

\details{

  \code{frscv} computes cross-validation for a regression spline
  estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous and ordinal/nominal factor
  predictors. The optimal \code{K}/\code{I} combination is returned
  along with other results (see below for return values).

  For the continuous predictors the regression spline model employs the
  tensor product B-spline basis matrix for a multivariate polynomial
  spline via the B-spline routines in the GNU Scientific Library and the
  \code{\link{tensor.prod.model.matrix}} function in the \pkg{mgcv}
  package.

  For the ordinal/nominal predictors the regression spline model uses
  indicator basis functions.

  Note that, since the approach is built on top of the regression spline
  for each unique ordinal/nominal combination, the splines must be
  estimable for every such unique combination. When this is not the case
  the user must combine data in a manner such that this is achieved.
  
}

\value{
  
  \code{frscv} returns a \code{crscv} object. Furthermore, the function
  \code{\link{summary}} supports objects of this type. The returned
  objects have the following components:
  
    \item{K}{ scalar/vector containing optimal degree(s) of spline }
    \item{I}{ scalar/vector containing an indicator of whether the
    predictor is included or not for each dimension of the
    ordinal/nominal factors }
    \item{K.mat}{ vector/matrix of values of \code{K} evaluated during search }  
    \item{basis.maxdim}{ maximum degree of spline }
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

## Compute the optimal K and I

cv <- frscv(x=xdata,y=y)
summary(cv)
}
\keyword{nonparametric}
