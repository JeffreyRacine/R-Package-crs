\name{plot.crs}
\alias{plot.crs}
\title{Plot Categorical Regression Spline Fits}
\description{
Plot fitted conditional mean, conditional quantile, gradient, and surface
displays for \code{\link{crs}} objects using the modern \pkg{np}-style plot
interface.
}
\usage{
\method{plot}{crs}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"crs"}.}
  \item{\dots}{plot controls. Common controls include
    \code{errors=c("none","asymptotic","bootstrap")}, \code{band},
    \code{alpha}, \code{bootstrap}, \code{B}, \code{data_overlay},
    \code{data_rug}, \code{legend}, \code{gradients},
    \code{gradient_order}, \code{perspective}, \code{renderer},
    \code{view}, \code{neval}, \code{output}, \code{np_grid_control},
    \code{np_boot_control}, and \code{np_render_control}.}
}
\details{
The default display is the fitted function, matching the modern
\pkg{np} regression plot convention. One-dimensional fits are shown as
partial fitted curves. For two continuous predictors,
\code{perspective=TRUE} requests a fitted surface, and
\code{renderer="rgl"} uses the \pkg{rgl} renderer when available.

Data overlays are controlled by \code{data_overlay}; support rugs are
controlled by \code{data_rug}. Intervals are requested with
\code{errors}; \code{band="all"} draws pointwise, simultaneous, and
Bonferroni bands and, by default, a legend. Set \code{legend=FALSE} to
suppress the legend or pass a legend position/list to customize it.

Use \code{output="data"} to return the plotting payload without drawing,
or \code{output="plot-data"} to draw and return the payload.
}
\value{
Invisibly returns the plotted data when \code{output="data"} or
\code{output="plot-data"} is requested; otherwise returns the plotted
object invisibly.
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
  \code{\link{crs}}, \code{\link{np_plot_controls}},
  \code{\link{np_grid_control}}, \code{\link{np_boot_control}},
  \code{\link{np_render_control}}
}
\examples{
set.seed(42)
n <- 100
x <- runif(n)
y <- sin(2*pi*x) + rnorm(n, sd = 0.25)
model <- crs(y ~ x, cv = "none", degree = 3, segments = 1)

plot(model)
plot(model, data_rug = TRUE)

\dontrun{
plot(model, errors = "bootstrap", band = "all", B = 99)
}
}
