\name{crs-package}
\alias{crs-package}
\docType{package}
\title{Categorical Regression Splines}
\description{

	This package provides a method for nonparametric regression that
	combines the (global) approximation power of regression splines for
	continuous predictors (\sQuote{\code{x}}) with the (local) power of
	kernel methods for categorical predictors (\sQuote{\code{z}}). The
	user also has the option of instead using indicator bases for the
	categorical predictors (default). When the predictors contain both
	continuous and categorical (discrete) data types, both approaches
	offer more efficient estimation than the traditional sample-splitting
	(i.e. \sQuote{frequency}) approach where the data is first broken into
	subsets governed by the categorical \code{z}.

  To cite the \pkg{crs} package type: \sQuote{\code{citation("crs")}}
  (without the single quotes).
  
  For a listing of all routines in the \pkg{crs} package type:
  \sQuote{\code{library(help="crs")}}.

  For a \sQuote{\code{\link{vignette}}} that presents an overview of the
  \pkg{crs} package type: \sQuote{\code{vignette("crs")}}.
  
}

\details{

  For the continuous predictors the regression spline model employs the
  B-spline basis matrix and for interactions may append the tensor
  product B-spline basis matrix for a multivariate polynomial spline via
  the B-spline routines in the GNU Scientific Library and the
  \code{\link{tensor.prod.model.matrix}} function. The B-spline function
  is the maximally differentiable interpolative basis function (B-spline
  stands for \dQuote{basis}-spline).
  
  For the discrete predictors the product kernel function is of the
  `Li-Racine' type (see Li and Racine (2007) for details) which is
  formed by constructing products of one of the following univariate
  kernels:

	\describe{

    \item{(\eqn{z} is discrete/nominal)}{
      \eqn{l(z_i,z,\lambda) = 1 }{l(z[i],z,lambda) =
        1} if \eqn{z_i=z}{z[i] = z}, and
      \eqn{\lambda}{lambda} if \eqn{z_i \neq z}{z[i] != z}. Note that
      \eqn{\lambda}{lambda} must lie between \eqn{0} and \eqn{1}.
      }

    \item{(\eqn{z} is discrete/ordinal)}{
      \eqn{l(z_i,z,\lambda) = 1}{l(z[i],z,lambda) = 1} if
      \eqn{|z_i-z|=0}{|z[i] - z| = 0}, and
      \eqn{\lambda^{|z_i-z|}}{lambda^|z_i-z|} if \eqn{|z_i -
        z|\ge1}{|z[i] - z|>=1}. Note that \eqn{\lambda}{lambda} must lie
      between \eqn{0} and \eqn{1}.
    }

  }

  Alternatively, for the ordinal/nominal predictors the regression
  spline model will use indicator basis functions.

}

\section{Usage Issues}{

  In order for \code{\link{plot}} to function properly with the options
  \code{mean=TRUE}, \code{deriv=TRUE}, and \code{ci=TRUE} you need to
  properly cast your categorical predictors prior to listing them in the
  \code{\link{formula}} interface.
  
}

\author{
  
  Jeffrey S. Racine <racinej@mcmaster.ca> and Zhenghua Nie <niez@mcmaster.ca>

  Maintainer: Jeffrey S. Racine <racinej@mcmaster.ca>

  I would like to gratefully acknowledge support from the Natural
  Sciences and Engineering Research Council of Canada
  (\url{http://www.nserc.ca}), the Social Sciences and Humanities
  Research Council of Canada (\url{http://www.sshrc.ca}), and the Shared
  Hierarchical Academic Research Computing Network
  (\url{http://www.sharcnet.ca}).
  
}
\references{

  Ma, S. and J.S. Racine and L. Yang (2011), \dQuote{Spline Regression
  in the Presence of Categorical Predictors,} manuscript.

  Ma, S. and J.S. Racine (2011), \dQuote{Additive Regression Splines
  With Irrelevant Categorical and Continuous Regressors,} manuscript.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
		Theory and Practice,} Princeton University Press.

}
\keyword{package}
