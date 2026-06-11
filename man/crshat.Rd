\name{crshat}
\alias{crshat}
\alias{crshat.crs}
\alias{print.crshat}
\title{CRS Hat Operator}
\usage{
crshat(object, ...)

\method{crshat}{crs}(object, newdata = NULL, y = NULL,
  output = c("matrix", "apply", "constraint"), deriv = 0,
  deriv.index = 1, rcond.min = 1e-08, use.svd.fallback = TRUE, ...)
}
\arguments{
  \item{object}{A fitted mean-regression object of class \code{"crs"}.}
  \item{newdata}{Optional evaluation data. If omitted, the training data are
    used.}
  \item{y}{Optional response vector or matrix. Required for
    \code{output="constraint"}; when supplied with \code{output="apply"}, the
    operator is applied directly without returning a dense hat matrix.}
  \item{output}{One of \code{"matrix"}, \code{"apply"}, or
    \code{"constraint"}.}
  \item{deriv, deriv.index}{Derivative controls reserved for future use.
    Currently only \code{deriv=0} is supported.}
  \item{rcond.min}{Minimum reciprocal condition number used by the weighted
    least-squares primitive before falling back to QR/SVD policy.}
  \item{use.svd.fallback}{Logical value controlling the weighted
    least-squares fallback policy.}
  \item{...}{Additional arguments, currently unused.}
}
\description{
Construct or apply the linear hat operator for fixed-structure CRS mean
regression fits. The helper is intended for fixed-\eqn{X} workflows such as
fast plot bootstrapping and constrained-estimation diagnostics.
}
\details{
For a fitted CRS mean-regression object, \code{crshat()} returns the linear
operator mapping the training response to fitted values at \code{newdata}. The
default \code{output="matrix"} returns a dense matrix with class
\code{c("crshat", "matrix")}. The \code{output="apply"} route applies the same
operator to \code{y} and is preferable when many right-hand sides are needed.

The helper preserves the fitted object's fixed spline/kernel structure,
including basis family, degree, segments, pruning state, categorical-kernel
bandwidths, weights, and the package's weighted least-squares rank policy. It
does not run bandwidth or degree selection.

Quantile CRS objects and derivative operators are intentionally unsupported at
present and fail with explicit errors.
}
\value{
For \code{output="matrix"}, a matrix of dimension
\code{nrow(newdata) x nrow(training data)}. For \code{output="apply"}, fitted
values produced by applying the operator to \code{y}. For
\code{output="constraint"}, the transpose operator scaled by the supplied
one-column response vector.
}
\seealso{
  \code{\link{crs}}, \code{\link{predict}}
}
\examples{
set.seed(42)
x <- runif(30)
y <- sin(2 * pi * x) + rnorm(30, sd = 0.1)
fit <- crs(y ~ x, cv = "none", degree = 2, segments = 1,
           display.warnings = FALSE, display.nomad.progress = FALSE)
nd <- data.frame(x = seq(0.1, 0.9, length.out = 5))
H <- crshat(fit, newdata = nd)
drop(H \%*\% fit$y)
crshat(fit, newdata = nd, output = "apply")
}
