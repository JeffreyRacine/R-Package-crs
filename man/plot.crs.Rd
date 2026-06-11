\name{plot.crs}
\alias{plot.crs}
\title{Plot Categorical Regression Spline Fits}
\description{
  Plot fitted conditional mean, conditional quantile, gradient, and
  bivariate surface displays for \code{\link{crs}} objects using the
  modern \pkg{np}-style plot interface.
}
\usage{
\method{plot}{crs}(x, \dots)
}
\arguments{
\subsection{Plot Object}{
This argument identifies the fitted CRS object to plot.
}

  \item{x}{
    an object of class \code{"crs"} returned by \code{\link{crs}}.
    Mean-regression objects plot fitted conditional means by default.
    Objects fit with \code{tau} plot fitted conditional quantiles.
  }

\subsection{Plot Options Passed Through \code{\dots}}{
Named options passed via \code{\dots} include:
}

    \strong{Data and evaluation controls}
    \describe{
      \item{\code{grid_control}}{
        optional detailed grid controls created by
        \code{\link{np_grid_control}}.
      }
      \item{\code{neval}}{
        the number of evaluation points used for each plotted continuous
        margin. For bivariate surfaces, the plotted grid has
        \code{neval * neval} evaluation points.
      }
      \item{\code{xq}}{
        the conditioning quantile used to hold non-plotted continuous
        variables fixed when drawing one-dimensional slices.
      }
      \item{\code{xtrim}}{
        trimming/expansion control for the evaluation support. The
        default evaluates over the observed support of the plotted
        continuous variable.
      }
    }

    \strong{Plot content controls}
    \describe{
      \item{\code{common_scale}}{
        logical; whether multiple one-dimensional panels should share a
        common vertical scale when feasible.
      }
      \item{\code{gradient_order}, \code{gradients}}{
        derivative plotting controls. \code{gradients=TRUE} requests
        gradients/effects rather than the fitted function, and
        \code{gradient_order} selects the derivative order. Current CRS
        gradient plotting requires a scalar \code{gradient_order}.
      }
      \item{\code{layout}}{
        \code{"auto"} or \code{TRUE} requests automatic
        \code{par(mfrow=...)} management for multi-panel plots.
        \code{"current"} or \code{FALSE} preserves the user's current
        graphics layout.
      }
      \item{\code{output}}{
        one of \code{"plot"}, \code{"data"}, or \code{"both"}.
        \code{"both"} is the public spelling for the historical
        \code{"plot-data"} behavior. \code{"data"} returns the plotting
        payload without drawing.
      }
    }

    \strong{Interval and bootstrap controls}
    \describe{
      \item{\code{alpha}}{
        the nominal significance level used for interval/band
        construction.
      }
      \item{\code{B}}{
        the bootstrap replication count. The default is \code{1999}.
      }
      \item{\code{band}}{
        one of \code{"pmzsd"}, \code{"pointwise"},
        \code{"bonferroni"}, \code{"simultaneous"}, or \code{"all"}.
        \code{band="all"} draws pointwise, simultaneous, and Bonferroni
        bands with a package-drawn legend unless suppressed.
      }
      \item{\code{boot_control}}{
        optional detailed bootstrap controls created by
        \code{\link{np_boot_control}}. This controls non-fixed
        bootstrap behavior, wild-bootstrap multipliers, and dependent
        bootstrap block length.
      }
      \item{\code{bootstrap}}{
        the bootstrap resampling scheme. Mean CRS plots support
        \code{"wild"}, \code{"inid"}, \code{"fixed"}, and
        \code{"geom"}. Quantile CRS plots default to \code{"inid"} and
        support \code{"inid"}, \code{"fixed"}, and \code{"geom"};
        \code{"wild"} is mean-only.
      }
      \item{\code{center}}{
        interval centering control. CRS currently supports
        \code{center="estimate"} only and fails fast for unsupported
        centering requests.
      }
      \item{\code{errors}}{
        one of \code{"none"}, \code{"bootstrap"}, or
        \code{"asymptotic"}.
      }
    }

    \strong{Rendering and view controls}
    \describe{
      \item{\code{data_overlay}}{
        logical; when \code{TRUE}, overlay observed response data on
        fitted mean/quantile displays where meaningful and expand plot
        limits to include the data. Overlay points default to small,
        transparent markers and can be customized with ordinary
        \code{points()} graphical arguments passed through \code{\dots}.
      }
      \item{\code{data_rug}}{
        logical; when \code{TRUE}, draw rug-style support cues for the
        plotted coordinates. For ordinary one-dimensional base plots this
        behaves like a standard rug. For eligible bivariate surface
        plots, the base renderer draws floor-rug ticks and
        \code{renderer="rgl"} draws analogous 3D floor-rug ticks.
      }
      \item{\code{legend}}{
        legend control for legends drawn by the package, including
        \code{band="all"} uncertainty legends. \code{TRUE} draws the
        default legend, \code{FALSE}, \code{NULL}, or \code{NA}
        suppresses it, a single character string supplies the legend
        position, and a list is merged into the default
        \code{\link{legend}} call for base plots or the default
        \code{rgl::legend3d()} call for \code{renderer="rgl"}. For
        routes with more than one package-drawn legend, named sublists
        such as \code{legend=list(bands=FALSE)} may be used when a route
        distinguishes those legend roles.
      }
      \item{\code{perspective}}{
        logical; whether two-continuous-predictor CRS fits should be
        displayed as bivariate fitted-function surfaces.
      }
      \item{\code{phi}, \code{theta}, \code{view}}{
        view controls for bivariate continuous surface plots. As in
        \pkg{np}, the \code{rgl} backend remaps the inherited base
        defaults \code{theta=0}, \code{phi=20} to a more informative
        oblique view.
      }
      \item{\code{render_control}}{
        optional detailed rendering controls created by
        \code{\link{np_render_control}}, including band-versus-bar
        uncertainty display controls for curve routes.
      }
      \item{\code{renderer}}{
        the surface-rendering backend for eligible bivariate continuous
        plots. \code{"base"} uses \code{\link{persp}}-style graphics.
        \code{"rgl"} uses the suggested \pkg{rgl} package and displays
        through the active \pkg{rgl} backend. On supported surface
        routes, \code{data_overlay=TRUE} adds observed data as an
        interactive point cloud, \code{data_rug=TRUE} adds 3D floor-rug
        ticks, and \code{errors!="none"} adds interval surfaces. With
        \code{band="all"}, CRS draws a matching interactive legend using
        \code{rgl::legend3d()}.
      }
    }

    \strong{Graphical parameters}
    \describe{
      \item{\code{border}, \code{cex}, \code{cex.axis}, \code{cex.lab},
      \code{cex.main}, \code{cex.sub}, \code{col}, \code{lty},
      \code{lwd}, \code{main}, \code{pch}, \code{sub}, \code{type},
      \code{xlab}, \code{xlim}, \code{ylab}, \code{ylim}, \code{zlab},
      \code{zlim}}{
        standard graphical controls forwarded to the underlying plotting
        primitives when relevant.
      }
      \item{\code{rgl.persp3d.<arg>}, \code{rgl.view3d.<arg>},
      \code{rgl.par3d.<arg>}, \code{rgl.grid3d.<arg>},
      \code{rgl.widget.<arg>}, \code{rgl.surface3d.<arg>},
      \code{rgl.legend3d.<arg>}}{
        backend-specific \pkg{rgl} controls. For example,
        \code{rgl.legend3d.bty="o"} restores a box around the
        \code{band="all"} \pkg{rgl} legend, and
        \code{rgl.surface3d.alpha} customizes interval-surface
        transparency.
      }
    }

    \strong{Reproducibility controls}
    \describe{
      \item{\code{random.seed}}{
        an optional random seed used to make bootstrap plotting
        reproducible.
      }
    }

\subsection{Additional Arguments}{
Further plotting controls are passed through \code{\dots} to the relevant
plot route.
}

  \item{\dots}{
    additional arguments supplied to control plotting behavior or passed
    through to underlying plotting helpers where supported. Supported
    named plotting controls are summarized above; standard graphical
    parameters are forwarded where relevant.
  }
}
\details{
  The default display is the fitted function, matching the modern
  \pkg{np} regression plot convention. One-dimensional fits are shown as
  partial fitted curves. For two continuous predictors,
  \code{perspective=TRUE} requests a fitted surface.

  Plot interval/error types are:
  \preformatted{
  "pmzsd"        : point estimate +/- z_(1-alpha/2) * standard error
  "pointwise"    : per-point two-sided interval
  "bonferroni"   : pointwise interval with alpha adjusted for m grid points
  "simultaneous" : joint-band interval from the bootstrap route
  }
  where \eqn{m} is the number of evaluation points used in the plotted
  curve or surface. \code{band="all"} overlays pointwise, simultaneous,
  and Bonferroni bands.

  For mean-regression plots, \code{errors="bootstrap"} defaults to the
  fast fixed-design wild bootstrap used by \pkg{np}, with
  \code{bootstrap="wild"} and \code{B=1999}. The refit bootstrap remains
  available with \code{bootstrap="inid"}. Mean plots also support the
  \pkg{np}-style block refit bootstraps \code{bootstrap="fixed"} and
  \code{bootstrap="geom"}; use \code{np_boot_control(blocklen=...)} to
  set the block length explicitly, otherwise the \pkg{np} automatic
  block-length selector is used.

  Quantile CRS plots, i.e. objects fit with \code{tau}, follow the
  \pkg{npqreg} convention and default to the refit
  \code{bootstrap="inid"} selector. \code{bootstrap="fixed"} and
  \code{bootstrap="geom"} are also available for quantile fitted-curve
  plots, while \code{bootstrap="wild"} remains mean-only.

  Unsupported legacy CRS plot arguments such as \code{ci},
  \code{deriv}, and \code{mean} now fail fast. New public plot controls
  use the \pkg{np}-style names documented here.
}
\value{
  Invisibly returns the plotted object for ordinary plotting calls.
  When \code{output="data"} or \code{output="both"} is requested, returns
  the plotting payload.
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
plot(model, errors = "bootstrap", band = "all")
plot(model, errors = "bootstrap", bootstrap = "inid", B = 99)
plot(model, errors = "bootstrap", bootstrap = "fixed", B = 99,
     boot_control = np_boot_control(blocklen = 5))

z <- runif(n)
y2 <- sin(2*pi*x) + z + rnorm(n, sd = 0.25)
smodel <- crs(y2 ~ x + z, cv = "none", degree = c(3, 1),
              segments = c(1, 1))
plot(smodel, perspective = TRUE, data_rug = TRUE)
plot(smodel, perspective = TRUE, renderer = "rgl",
     errors = "bootstrap", band = "all", B = 99)

qmodel <- crs(y ~ x, cv = "none", degree = 3, segments = 1, tau = 0.5)
plot(qmodel, errors = "bootstrap", B = 99)
}
}
