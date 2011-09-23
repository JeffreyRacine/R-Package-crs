gsl.bs <- function(...) UseMethod("gsl.bs")

gsl.bs.default <- function(x,
                           degree=3,
                           nbreak=2,
                           deriv=0,
                           x.min=NULL,
                           x.max=NULL,
                           intercept=FALSE,
                           knots=NULL,
                           ...) {

  x <- as.vector(x)
  ## Some error checking

  if(degree <= 0) stop(" degree must be a positive integer")
  if(deriv < 0) stop(" deriv must be a non-negative integer")
  if(deriv > degree +1 ) stop(" deriv must be smaller than degree plus 2") 
  if(nbreak <= 1) stop(" nbreak must be at least 2")

#  if(!is.null(knots)) nbreak <- length(knots)
  if(!is.null(knots)&&length(knots)!=nbreak) {
    nbreak <- length(knots)
    warning(paste(" nbreak and knots vector do not agree: resetting nbreak to", nbreak))
  }

  ## For evaluation (newx) must use min/max for x unless otherwise
  ## specified - check that mix < max


  if(!is.null(x.min)&!is.null(x.max)) if(x.min >= x.max) stop(" x.min must be less than x.max")
  if(is.null(x.min)) {
			x.min <- min(x)
			ol <- FALSE
	}
	else {
			ol <- x < x.min 
	}
  if(is.null(x.max)) {
			x.max <- max(x)
			or <- FALSE
	}
	else {
			or <- x > x.max 
	}

  ## we should check Xs  are within knots   
## It seems that 'knots' is given,  we don't require x.min and x.max,  otherwise,  we check the boundary.

	if(!is.null(knots)) {
			ol <- ol | x < min(knots)
			or <- or | x > max(knots)
	}

	outside <- ol | or 

	if(any(outside)){
			warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
			ord <- degree + 1
#			Boundary.knots <- c(x.min, x.max)
#			Aknots <- sort(c(rep(Boundary.knots, degree+1), knots))
#			nbreak1 <- length(Aknots) 
			derivs<- 0:degree
			scalef <- gamma(1L:ord)
			B <- array(0, c(length(x), nbreak+degree-1))
			
			if(any(ol)){##the following code follows bs,  but it does not provide the derivatives,  we need to know how to compute derivatives when
					# these happen
					k.pivot <- x.min
					xl <- cbind(1, outer(x[ol]-k.pivot, 1L:degree,  "^"))
					tt <- bs.des(rep(k.pivot, ord), degree, nbreak, deriv=derivs, x.min, x.max)  #If 'knots' is given,  we have not implemented it.
					B[ol, ] <- xl %*% (tt/scalef)
			}
			if(any(or)){
					k.pivot <- x.max
					xr <- cbind(1, outer(x[or]-k.pivot, 1L:degree,  "^"))
					tt <- bs.des(rep(k.pivot, ord), degree, nbreak,deriv= derivs, x.min, x.max)
					B[or, ] <- xr %*% (tt/scalef)
			}
			if(any(inside <- !outside))
					B[inside, ] <- bs.des(x[inside], degree, nbreak, rep(deriv, length(x[inside])), x.min, x.max, knots)
	}
	else {
			B <- bs.des(x, degree, nbreak, rep(deriv, length(x)), x.min, x.max,knots)
	}


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

bs.des     <- function(x,
											 degree=3,
											 nbreak=2,
											 deriv=integer(length(x)),   ## now deriv is a vector to represent the derivivative for each element of x 
											 x.min=NULL,
											 x.max=NULL,
											 knots=NULL,
											 ...) {

		x <- as.vector(x)
		n <- length(x)
		deriv <- as.vector(deriv)

		if(!is.null(x.min)&!is.null(x.max)) if(x.min >= x.max) stop(" x.min must be less than x.max")
		if(is.null(x.min)) x.min <- min(x)
		if(is.null(x.max)) x.max <- max(x)

		## 0 == don't use user supplied knots, 1 = use

		knots.int <- ifelse(is.null(knots), 0, 1)

		ncol <- nbreak+degree-1;

		if(all(deriv==0)) {

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
										as.integer(max(deriv)), #for setting the size of memory
										as.double(x.min),
										as.double(x.max),
										as.double(knots),                
										as.integer(knots.int),
										Bx = double(n*ncol),
										PACKAGE="crs" )

		}

		B <- matrix(data=myout$Bx, nrow = n, ncol = ncol, byrow = TRUE)

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

#				if(min(newx)<x.min || max(newx)>x.max) {
#						warning(" evaluation data lies beyond spline support: resetting those values to min/max")
#						newx[newx < x.min] <- x.min
#						newx[newx > x.max] <- x.max
#						newx.ind <- sort(c(which(newx < x.min),which(newx > x.max)))
#				}

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
		attr(B, "newx.trimmed") <- newx.ind

		return(B)

}
