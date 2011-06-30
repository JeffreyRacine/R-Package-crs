\name{npglpreg}
\alias{npglpreg}
\alias{npglpreg.default}
\alias{npglpreg.formula}

\title{Generalized Local Polynomial Regression}
\description{
  
  \code{npglpreg} computes a generalized local polynomial kernel
    regression estimate of a one (1) dimensional dependent variable on
    an \code{r}-dimensional vector of continuous and categorical
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
         \dots)

\method{npglpreg}{formula}(formula,
         data = list(),
         tydat = NULL,
         txdat = NULL,
         eydat = NULL,
         exdat = NULL,
         bws = NULL,
         degree = NULL,
         degree.max = 5,
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
         \dots)


}

\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{tydat}{
%%     ~~Describe \code{tydat} here~~
}
  \item{txdat}{
%%     ~~Describe \code{txdat} here~~
}
  \item{eydat}{
%%     ~~Describe \code{eydat} here~~
}
  \item{exdat}{
%%     ~~Describe \code{exdat} here~~
}
  \item{bws}{
%%     ~~Describe \code{bws} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{degree.max}{
%%     ~~Describe \code{degree.max} here~~
}
  \item{leave.one.out}{
%%     ~~Describe \code{leave.one.out} here~~
}
  \item{ukertype}{
%%     ~~Describe \code{ukertype} here~~
}
  \item{okertype}{
%%     ~~Describe \code{okertype} here~~
}
  \item{bwtype}{
%%     ~~Describe \code{bwtype} here~~
}
  \item{cv}{
%%     ~~Describe \code{cv} here~~
}
  \item{cv.func}{
%%     ~~Describe \code{cv.func} here~~
}
  \item{opts}{
%%     ~~Describe \code{opts} here~~
}
  \item{nmulti}{
%%     ~~Describe \code{nmulti} here~~
}
  \item{random.seed}{
%%     ~~Describe \code{random.seed} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
  This function is in beta status until further notice.

  Important note - prediction with newdata only works if there exists
  more than one row in the newdata data frame.
}
\value{
TBA
}
\references{
TBA
}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
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

model.glp <- npglpreg(y~x1+x2)
}
\keyword{ nonparametric }
\keyword{ regression }

