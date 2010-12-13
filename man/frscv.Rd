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
      basis.maxdim = 5,
      complexity = c("degree-knots","degree","knots"),
      knots = c("quantiles","uniform"),
      basis = c("additive-tensor","additive","tensor","auto"),
      cv.norm = c("L2","L1"),
      degree = degree,
      segments = segments)
}

\arguments{

  \item{y}{
    continuous univariate vector
  }
  
  \item{xz}{
    continuous and/or ordinal/nominal predictors
  }

  \item{basis.maxdim}{ the maximum dimension of the B-spline basis for
    each of the continuous predictors (default \code{basis.maxdim=5})}

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
  equally spaced intervals. }

\item{basis}{ a character string (default
    \code{basis="additive-tensor"}) indicating whether the
    additive/tensor, additive, or combined tensor product B-spline basis
    matrix for a multivariate polynomial spline should be used. Note
    this can be automatically determined by cross-validation if
    \code{cv=TRUE} and \code{basis="auto"}, and is an \sQuote{all or
    none} proposition (i.e. interaction terms for all predictors or for
    no predictors given the nature of \sQuote{tensor products}). Note
    also that if there is only one predictor this defaults to
    \code{basis="additive"} to avoid unnecessary computation as the
    spline bases are equivalent in this case }

  \item{cv.norm}{ a character string (default \code{cv.norm="L2"}
  indicating whether the L2-norm is used (sum of squared delete-one
  residuals) or the L1-norm (sum of absolute delete-one-residuals) in
  the cross-validation function }

  \item{degree}{ integer/vector specifying the degree of the B-spline
  basis for each dimension of the continuous \code{x}}
  
  \item{segments}{ integer/vector specifying the number of segments of
  the B-spline basis for each dimension of the continuous \code{x}
  (i.e. number of knots minus one)}

}

\details{

  \code{frscv} computes cross-validation for a regression spline
  estimate of a one (1) dimensional dependent variable on an
  \code{r}-dimensional vector of continuous and ordinal/nominal factor
  predictors. The optimal \code{K}/\code{I} combination is returned
  along with other results (see below for return values).

  For the continuous predictors the regression spline model employs
  either the additive, tensor product, or combined additive/tensor
  product B-spline basis matrix for a multivariate polynomial spline via
  the B-spline routines in the GNU Scientific Library and the
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
  
  \item{K}{ scalar/vector containing optimal degree(s) of spline or
    number of segments }
  \item{I}{ scalar/vector containing an indicator of whether the
    predictor is included or not for each dimension of the
    ordinal/nominal factors }
  \item{K.mat}{ vector/matrix of values of \code{K} evaluated during search }  
  \item{basis.maxdim}{ maximum degree of B-spline basis }
  \item{cv.func}{ objective function value at optimum }
  \item{cv.func.vec}{ vector of objective function values at each degree
    of spline or number of segments in \code{K.mat}}

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

## Compute the optimal K and I, determine optimal number of knots, set
## spline degree for x to 3

cv <- frscv(x=xdata,y=y,complexity="knots",degree=c(3))
summary(cv)
}
\keyword{nonparametric}
