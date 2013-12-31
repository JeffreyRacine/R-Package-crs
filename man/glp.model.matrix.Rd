\name{glp.model.matrix}
\alias{glp.model.matrix}
\title{Utility function for constructing generalized polynomial smooths}
\description{
  Produce model matrices for a generalized polynomial smooth from the
  model matrices for the marginal bases of the smooth.
}
\usage{
glp.model.matrix(X)
}
\arguments{
  \item{X}{a list of model matrices for the marginal bases of a smooth} 
}
\details{This function computes a generalized polynomial where the
  orders of each term entering the polynomial may vary.}

\value{A model matrix for a generalized polynomial smooth.}

\references{

  Hall, P. and J.S. Racine (forthcoming), \dQuote{Cross-Validated
  Generalized Local Polynomial Regression,} Journal of Econometrics.

}
\author{ Jeffrey S. Racine \email{racinej@mcmaster.ca}} 

\examples{
X <- list(matrix(1:4,2,2),matrix(5:10,2,3))
glp.model.matrix(X)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..




