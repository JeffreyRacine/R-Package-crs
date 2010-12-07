gsl.bs <- function(...) UseMethod("gsl.bs")

gsl.bs.default <- function(x,
                           degree=3,
                           nbreak=2,
                           deriv=0,
                           x.min=NULL,
                           x.max=NULL,
                           intercept=FALSE,
                           knots=c("quantiles","uniform"),
                           ...) {

  x <- as.vector(x)
  n <- length(x)

  knots <- match.arg(knots)

  if(degree <= 0) stop(" degree must be a positive integer")
  if(deriv < 0) stop(" deriv must be a non-negative integer")
  if(nbreak <= 1) stop(" nbreak must be at least 2")

  ## For evaluation (newx) must use min/max for x unless otherwise
  ## specified

  if(is.null(x.min)) x.min <- min(x)
  if(is.null(x.max)) x.max <- max(x)

  ## Construct quantiles in case they are requested (trivial, no
  ## overhead so may as well compute rather than test etc.)
  
  q.vec <- quantile(x,seq(0,1,length=nbreak))

  ## Replace endpoints with x.min and x.max if overloaded, otherwise
  ## this has no effect (0 and 1 quantiles are sample min/max)
  
  q.vec[1] <- x.min
  q.vec[nbreak] <- x.max

  ## Integer to pass to C code, 0 uniform, 1 quantile

  knots.int <- ifelse(knots == "uniform", 0, 1)

  ncol <- nbreak+degree-1;

  if(deriv==0) {

    myout <- .C("gsl_bspline",
                as.double(x),
                as.integer(n),
                as.integer(degree),
                as.integer(nbreak),
                as.double(x.min),
                as.double(x.max),
                as.double(q.vec),
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
                as.double(q.vec),                
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
  attr(B, "class") <- c("gsl.bs","matrix")

  return(B)

}

predict.gsl.bs <- function(object,
                           newx=NULL,
                           ...) {

  newx.ind <- NULL

  if(is.null(newx))

    ## If No new data provided, return sample fit.
    B <- object

  else{

    x.min <- attr(object, "x.min")
    x.max <- attr(object, "x.max")

    newx <- as.numeric(newx)

    if(min(newx)<x.min || max(newx)>x.max) {
      warning(" evaluation data lies beyond spline support: resetting those values to min/max")
      newx[newx < x.min] <- x.min
      newx[newx > x.max] <- x.max
      newx.ind <- sort(c(which(newx < x.min),which(newx > x.max)))
    }

    B <- gsl.bs(newx,
                degree=attr(object, "degree"),
                nbreak=attr(object, "nbreak"),
                deriv=attr(object, "deriv"),
                intercept=attr(object, "intercept"),
                x.min=x.min,
                x.max=x.max)

  }

  attr(B, "newx") <- newx
  attr(B, "newx.trimmed") <- newx.ind

  return(B)

}
