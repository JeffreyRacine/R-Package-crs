\name{crs-package}
\alias{crs-package}
\docType{package}
\title{Nonparametric Regression Splines with Continuous and Categorical Predictors}
\description{

	This package provides a method for nonparametric regression that
	combines the (global) approximation power of regression splines for
	continuous predictors (\sQuote{\code{x}}) with the (local) power of
	kernel methods for categorical predictors (\sQuote{\code{z}}). The
	user also has the option of instead using indicator bases for the
	categorical predictors. When the predictors contain both continuous
	and categorical (discrete) data types, both approaches offer more
	efficient estimation than the traditional sample-splitting
	(i.e. \sQuote{frequency}) approach where the data is first broken into
	subsets governed by the categorical \code{z}.

  To cite the \pkg{crs} package type: \sQuote{\code{citation("crs")}}
  (without the single quotes).
  
  For a listing of all routines in the \pkg{crs} package type:
  \sQuote{\code{library(help="crs")}}.

  For a listing of all demos in the \pkg{crs} package type:
  \sQuote{\code{demo(package="crs")}}.

  For a \sQuote{\code{\link{vignette}}} that presents an overview of the
  \pkg{crs} package type: \sQuote{\code{vignette("crs")}}.
  
}

\details{

  For the continuous predictors the regression spline model employs the
  B-spline basis matrix using the B-spline routines in the GNU
  Scientific Library (\url{https://www.gnu.org/software/gsl/}).

   The \code{\link{tensor.prod.model.matrix}} function is used to
  construct multivariate tensor spline bases when \code{basis="tensor"}
  and uses additive B-splines otherwise (i.e. when
  \code{basis="additive"}).

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

\author{
  
  Jeffrey S. Racine \email{racinej@mcmaster.ca} and Zhenghua Nie \email{niez@mcmaster.ca}

  Maintainer: Jeffrey S. Racine \email{racinej@mcmaster.ca}

  I would like to gratefully acknowledge support from the Natural
  Sciences and Engineering Research Council of Canada
  (\url{https://www.nserc-crsng.gc.ca}), the Social Sciences and Humanities
  Research Council of Canada (\url{https://www.sshrc-crsh.gc.ca}), and the Shared
  Hierarchical Academic Research Computing Network
  (\url{https://www.sharcnet.ca}).
  
}
\references{

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
		Theory and Practice,} Princeton University Press.

  Ma, S. and J.S. Racine and L. Yang (2015), \dQuote{Spline
  Regression in the Presence of Categorical Predictors,} Journal of
  Applied Econometrics, Volume 30, 705-717.

  Ma, S. and J.S. Racine (2013), \dQuote{Additive Regression
  Splines with Irrelevant Categorical and Continuous Regressors,}
  Statistica Sinica, Volume 23, 515-541.

}
\keyword{package}
