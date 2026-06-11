\name{np_plot_controls}
\alias{np_plot_controls}
\alias{np_boot_control}
\alias{np_grid_control}
\alias{np_render_control}
\alias{crs_boot_control}
\alias{crs_grid_control}
\alias{crs_render_control}
\title{Plot Control Helpers}

\description{
Control objects for modern \pkg{crs} plot methods. The \code{np_*}
helper names mirror the plot interface used by \pkg{np}; the \code{crs_*}
names are package-local aliases that return compatible objects.
}

\usage{
np_boot_control(nonfixed = c("exact", "frozen"),
                wild = c("rademacher", "mammen"),
                blocklen = NULL)

np_grid_control(xtrim = NULL, xq = NULL, slices = NULL)

np_render_control(style = c("band", "bar"),
                  bar = c("|", "I"),
                  bar_num = NULL)

crs_boot_control(...)

crs_grid_control(...)

crs_render_control(...)
}

\arguments{
  \item{nonfixed}{Handling for nonfixed regressors in bootstrap plots.}
  \item{wild}{Wild-bootstrap multiplier family.}
  \item{blocklen}{Optional positive block length for
  \code{bootstrap="fixed"} and \code{bootstrap="geom"}. If omitted for
  mean \code{\link{crs}} plots, CRS uses the same automatic
  block-length selector as \pkg{np}.}
  \item{xtrim}{Optional length-two trimming quantiles for plot grids.}
  \item{xq}{Optional conditioning quantile or vector of conditioning quantiles.}
  \item{slices}{Reserved for future explicit slice selection.}
  \item{style}{Interval rendering style.}
  \item{bar}{Error-bar glyph style.}
  \item{bar_num}{Optional number of error bars to draw.}
  \item{\dots}{Arguments forwarded to the corresponding \code{np_*} helper.}
}

\value{
A small S3 list consumed by \code{plot.crs()} and related plot methods.
}

\seealso{
\code{\link{plot.crs}}, \code{\link{crs}}
}

\examples{
boot <- np_boot_control(nonfixed = "frozen", wild = "mammen")
grid <- np_grid_control(xtrim = c(0.05, 0.95), xq = 0.5)
render <- np_render_control(style = "bar", bar = "I")
}
