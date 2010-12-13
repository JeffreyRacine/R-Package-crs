gsl.bs <- function(...) UseMethod("gsl.bs")

gsl.bs.default <- function(x,
                           degree=3,
                           nbreak=NULL,
                           deriv=0,
                           x.min=NULL,
                           x.max=NULL,
                           intercept=FALSE,
                           knots=NULL,
                           ...) {

  x <- as.vector(x)
  n <- length(x)

  ## Some error checking

  if(degree <= 0) stop(" degree must be a positive integer")
  if(deriv < 0) stop(" deriv must be a non-negative integer")
  if(!is.null(nbreak)&&(nbreak <= 1)) stop(" nbreak must be at least 2")
  if(is.null(knots)&&is.null(nbreak)) stop(" either knots or nbreak must be provided")

  ## We may want to disable warnings as these things might be standard?

  if(!is.null(knots)&&is.null(nbreak)) nbreak <- length(knots)
  if(!is.null(knots)&!is.null(nbreak)&&length(knots)!=nbreak) {
    nbreak <- length(knots)
    warning(paste(" nbreak and knots vector do not agree: resetting nbreak to", nbreak))
  }

  ## For evaluation (newx) must use min/max for x unless otherwise
  ## specified - check that mix < max

  if(!is.null(x.min)&!is.null(x.max)) if(x.min >= x.max) stop(" x.min must be less than x.max")
  if(is.null(x.min)) x.min <- min(x)
  if(is.null(x.max)) x.max <- max(x)

  ## This appears to replicate the behaviour of bs() in the splines
  ## package when knots are supplied but x lies outside of the knot
  ## intervals.

  if(!is.null(knots)) {
    if(min(x) < min(knots)) {
      knots <- c(x.min,knots)
      nbreak <- length(knots)
      warning(" x.min < min(knots): extending knot range for out-of-support evaluation")    
    }
    if(x.max > max(knots)) {
      knots <- c(knots,x.max)
      nbreak <- length(knots)
      warning(" max(x) > max(knots): extending knot range for out-of-support evaluation")
    }
  }

  if(is.null(knots)) {
    if(x.min > min(x)) {
      x.min <- min(x)
      warning(" x.min > min(x): extending knot range for out-of-support evaluation")    
    }
    if(x.max < max(x)) {
      x.max <- max(x)
      warning(" x.max < max(x): extending knot range for out-of-support evaluation")
    }
  }

  ## 0 == don't use user supplied knots, 1 = use

  knots.int <- ifelse(is.null(knots), 0, 1)

  if(!is.null(knots)) {
    knots <- unique(sort(knots)) ## unique not sufficient?
    nbreak.unique <- length(knots)
    if(nbreak.unique < nbreak) {
      warning(" nbreak dynamically reduced due to non-uniqueness of quantile knot vector")
      nbreak <- nbreak.unique
    }
    if(nbreak <= 1) stop(" dynamically adjusted nbreak for quantile knot vector must be at least 2")
  }

  ncol <- nbreak+degree-1;

  if(deriv==0) {

    myout <- .C("gsl_bspline",
                as.double(x),
                as.integer(n),
                as.integer(degree),
                as.integer(nbreak),
                as.double(x.min),
                as.double(x.max),
                as.double(knots),
                as.integer(knots.int),
                Bx = double(n*ncol),
                PACKAGE="crs" )

  } else {

    myout <- .C("gsl_bspline_deriv",
                as.double(x),
                as.integer(n),
                as.integer(degree),
                as.integer(nbreak),
                as.integer(deriv),
                as.double(x.min),
                as.double(x.max),
                as.double(knots),                
                as.integer(knots.int),
                Bx = double(n*ncol),
                PACKAGE="crs" )

  }

  B <- matrix(data=myout$Bx, nrow = n, ncol = ncol, byrow = TRUE)
  if(!intercept) B <- B[,-1,drop=FALSE]

  attr(B, "degree") <- degree
  attr(B, "nbreak") <- nbreak
  attr(B, "deriv") <- deriv
  attr(B, "x.min") <- x.min
  attr(B, "x.max") <- x.max
  attr(B, "intercept") <- intercept
  attr(B, "knots") <- knots
  attr(B, "class") <- c("gsl.bs","matrix")

  return(B)

}

predict.gsl.bs <- function(object,
                           newx=NULL,
                           ...) {

  if(is.null(newx)) {

    ## If No new data provided, return sample fit.
    B <- object

  } else {

    x.min <- attr(object, "x.min")
    x.max <- attr(object, "x.max")

    newx <- as.numeric(newx)

    B <- gsl.bs(newx,
                degree=attr(object, "degree"),
                nbreak=attr(object, "nbreak"),
                deriv=attr(object, "deriv"),
                intercept=attr(object, "intercept"),
                knots=attr(object, "knots"),
                x.min=x.min,
                x.max=x.max)

  }

  attr(B, "newx") <- newx

  return(B)

}
