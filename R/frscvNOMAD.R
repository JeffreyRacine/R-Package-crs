# This function conducts factor regression spline
## cross-validation. It takes as input a data.frame xz containing a
## mix of numeric and factor predictors and a vector y. A range of
## arguments can be provided, and one can do search on both the degree
## and knots ("degree-knots") or the degree holding the number of
## knots (segments+1) constant or the number of knots (segments+1)
## holding the degree constant. Two basis types are supported
## ("additive", or "tensor") and the argument "auto" will choose the
## basis type automatically.


frscvNOMAD <- function(xz,
											 y,
											 degree.max=10, 
											 segments.max=10, 
											 degree.min=0, 
											 segments.min=1, 
											 complexity=c("degree-knots","degree","knots"),
											 knots=c("quantiles","uniform"),
											 basis=c("additive","tensor","auto"),
											 cv.func=c("cv.ls","cv.gcv","cv.aic"),
											 degree=degree,
											 segments=segments, 
											 include=include, 
											 opts=list(),
											 nmulti=0) {

		complexity <- match.arg(complexity)
		knots <- match.arg(knots)
		basis <- match.arg(basis)
		cv.func <- match.arg(cv.func)

		if ( missing(include) || is.null(include)) {
				include <- NULL
		}
		if(degree.min < 0 ) degree.min <- 0
		if(segments.min < 1 ) segments.min <- 1
		if(degree.max < degree.min) degree.max <- (degree.min + 1)
		if(segments.max < segments.min) segments.max <- (segments.min + 1)
		if(missing(degree)) degree <- NULL
		if(missing(segments)) segments <- NULL


		t1 <- Sys.time()

		cv.nomad <- function(x,
												 z=NULL,
												 y,
												 degree.max=degree.max, 
												 segments.max=segments.max, 
												 degree.min=degree.min, 
												 segments.min=segments.min, 
												 knots=knots,
												 basis=basis,
												 cv.func=cv.func, 
												 complexity=complexity, 
												 degree=degree, 
												 segments=segments,
												 include=include, 
												 opts=opts,
												 print.output=print.output, 
												 nmulti=nmulti) {

				if( missing(x) || missing(y) ) stop(" you must provide input, x, y")

				## Presumes x (continuous predictors) exist, but z
				## (ordinal/nominal factors) can be optional

				n <- length(y)
				num.x <- NCOL(x)

				if(NROW(y) != NROW(x)) stop(" x and y have differing numbers of observations")

				## K is a matrix, column 1 degree, column 2 segments, either or
				## both can be determined via cv so need to take care to allow
				## user to select knots (degree fixed), degree (knots fixed), or
				## both degree and knots. The values used to evaluate the cv
				## function are passed below.

				eval.cv <- function(input, params){

						complexity <- params$complexity
						segments <- params$segments
						degree <- params$degree
						include <- params$include
						x <- params$x
						y <- params$y
						z <- params$z
						knots <- params$knots
						cv.func <- params$cv.func
						basis <- params$basis

						num.x <- NCOL(x)
						num.z <- 0
						I <- include

						if(complexity=="degree-knots") {
								K <- round(cbind(input[1:num.x],input[(num.x+1):(2*num.x)]))
								if(!is.null(z)  ) {
										num.z <- NCOL(z)
										I <- round(input[(2*num.x+1):(2*num.x+num.z)])
								} 
						}	
						else if(complexity=="degree") {
								K<-round(cbind(input[1:num.x],segments))
								if(!is.null(z)  ) {
										num.z <- NCOL(z)
										I <- round(input[(num.x+1):(num.x+num.z)])
								} 
						}	
						else if(complexity=="knots")
						{
								K<-round(cbind(degree, input[1:num.x]))
								if(!is.null(z)  ) {
										num.z <- NCOL(z)
										I <- round(input[(num.x+1):(num.x+num.z)])
								} 
						}

						basis.opt <-  basis;
						if(basis=="auto"){
								basis.opt <-"additive"
								cv <- cv.factor.spline(x=x,
																			 y=y,
																			 z=z,
																			 K=K,
																			 I=I,
																			 knots=knots,
																			 basis=basis.opt,
																			 cv.func=cv.func)

								cv.tensor <- cv.factor.spline(x=x,
																							y=y,
																							z=z,
																							K=K,
																							I=I,
																							knots=knots,
																							basis="tensor",
																							cv.func=cv.func)
								if(cv>cv.tensor){
										cv <- cv.tensor
										basis.opt <-"tensor"
								}
						}
						else {
								cv <- cv.factor.spline(x=x,
																			 y=y,
																			 z=z,
																			 K=K,
																			 I=I,
																			 knots=knots,
																			 basis=basis.opt,
																			 cv.func=cv.func)
						}
						attr(cv, "basis.opt")<-basis.opt

						return(cv)
				}

				##generate params
				params <- list()
				params$complexity <- complexity
				params$segments <- segments
				params$degree <- degree
				params$include <- include
				params$x <- x
				params$y <- y
				params$z <- z
				params$knots <- knots
				params$cv.func <- cv.func
				params$basis <- basis

				# initial value
				if(!is.null(z)) {
						num.z <- NCOL(z)
				}
				else {
						num.z <-0
						include <- NULL
				}

				xsegments <- segments
				xdegree <- degree
				xinclude <- include

				if(is.null(xdegree)) xdegree <- sample(degree.min:degree.max, num.x, replace=T)
				if(is.null(xsegments)) xsegments <- sample(segments.min:segments.max, num.x, replace=T)
				if(is.null(xinclude)) xinclude <- sample(0:1, num.z, replace=T)
				

				if(complexity =="degree-knots") {
						x0 <- c(xdegree, xsegments, xinclude)
						bbin <-rep(1, num.x*2+num.z)
						#bounds segments cannot be samller than 2
						lb <- c(rep(degree.min,num.x), rep(segments.min, num.x), rep(0, num.z) )
						ub <- c(rep(degree.max, num.x), rep(segments.max,num.x), rep(1, num.z))
				}	
				else if(complexity=="degree") {
						x0 <- c(xdegree, xinclude)
						bbin <-rep(1, num.x+num.z)
						lb <- c(rep(degree.min, num.x),  rep(0, num.z) )
						ub <- c(rep(degree.max, num.x), rep(1, num.z))
				}	
				else if(complexity=="knots")
				{
						x0 <- c(xsegments, xinclude)
						bbin <-rep(1, num.x+num.z)
						#bounds segments cannot be samller than 2
						lb <- c(rep(segments.min, num.x),  rep(0, num.z) )
						ub <- c(rep(segments.max, num.x), rep(1, num.z))
				}

				if(length(x0) != length(lb)) stop(" x0 and bounds have differing numbers of variables")

				#no constraints
				bbout <-c(0)

				## This could be a passable parameter or set appropriately

				## Manual says preceed by r means relative to up and lb... not
				## quite what I was looking for

				solution<-snomadr(eval.f=eval.cv,
													n=length(x0),
													x0=as.numeric(x0),
													bbin=bbin,
													bbout=bbout,
													lb=lb,
													ub=ub,
													nmulti=as.integer(nmulti),
													opts=opts,
													print.output=print.output, 
													params=params);

				if(basis == "auto") 
				attr(solution, "basis.opt") <- attributes(eval.cv(solution$solution, params))$basis.opt

				return(solution)
		}

		console <- newLineConsole()
		print.output <- FALSE
		if(!is.null(opts$DISPLAY_DEGREE)){
				if(opts$DISPLAY_DEGREE>0){
						print.output <- TRUE
						console <- printPush("Being Solved by NOMAD...\n",console = console)
				}
		}
		else {
				print.output <- TRUE
				console <- printPush("Being Solved by NOMAD...\n",console = console)
		}
		## Take data frame x and parse into factors (z) and numeric (x)

		if(!is.data.frame(xz)) stop(" xz must be a data frame")

		xztmp <- splitFrame(xz)
		x <- xztmp$x
		z <- xztmp$z
		if(is.null(z)) {
				include <- NULL
				num.z <- 0
		} else {
				num.z <- NCOL(z)
				if(!is.null(include) && length(include)!= num.z){ 
						warning(paste(" the length of include (", length(include),") is not the same as the length of z (", num.z, ")",sep=""))
						include <- NULL
				}
				if (any(include <0) || any(include > 1)) include <- NULL
		}

		num.x <- ncol(x)
		n <- nrow(x)

		if(missing(x) || missing(y)) stop (" you must provide x and y")

		if(complexity=="degree") {
				if(missing(segments)||is.null(segments)) stop("segments missing for cross-validation of spline degree")
				if(length(segments)!=num.x) stop(" segments vector must be the same length as x")  
				if(!is.null(degree) && length(degree) == num.x) { #check initial values should be in the bounds
						if(any(degree < degree.min)||any(degree>degree.max)) {
								warning(paste(" The provided initial values for the degree are not in the bounds.", sep=""))
								degree <- NULL
						}
				}
				else
						degree <- NULL
		} else if(complexity=="knots") {
				if(missing(degree)||is.null(degree)) stop("degree missing for cross-validation of number of spline knots")
				if(length(degree)!=num.x) stop(" degree vector must be the same length as x")
				if(!is.null(segments) && length(segments) == num.x) {
						if(any(segments < segments.min)||any(segments > segments.max)) {
								warning(paste(" The provided initial values for the semgents are not in the bounds.", sep=""))
								segments <- NULL
						}
				}
				else
						segments <- NULL
		}
		else {
				if(!is.null(degree) && length(degree) == num.x) { #check initial values should be in the bounds
						if(any(degree < degree.min)||any(degree>degree.max)) {
								warning(paste(" The provided initial values for the degree are not in the bounds.", sep=""))
								degree <- NULL
						}
				}
				else {
						degree <- NULL
				}
				if(!is.null(segments) && length(segments) == num.x) {
						if(any(segments < segments.min)||any(segments > segments.max)) {
								warning(paste(" The provided initial values for the semgents are not in the bounds.", sep=""))
								segments <- NULL
						}
				}
				else
						segments <- NULL
		}

		## For factor regression spline, if there is only one predictor
		## (i.e. num.x + num.z = 1) disable auto, set to additive (tensor in
		## this case, so don't waste time doing both).

		if(num.x+num.z==1 & basis == "auto") basis <- "additive"

		if(degree.max < 1 || segments.max < 1 ) stop(" degree.max or segments.max must be greater than or equal to 1")

		## sovle by NOMAD
		nomad.solution <- cv.nomad(x,
															 z,
															 y,
															 degree.max=degree.max, 
															 segments.max=segments.max, 
															 degree.min=degree.min, 
															 segments.min=segments.min, 
															 knots=knots,
															 basis=basis,
															 cv.func=cv.func, 
															 complexity=complexity, 
															 degree=degree, 
															 segments=segments, 
															 include=include, 
															 opts=opts,
															 print.output=print.output, 
															 nmulti=nmulti) 

		# in crs,  we do not need the file best_x.txt
		if(nmulti > 0) file.remove("best_x.txt")

		t2 <- Sys.time()

		##output

		I.opt <- NULL 
		cv.min <- nomad.solution$objective
		basis.opt <- basis
		if(complexity=="degree-knots") {
				K.opt <- as.integer(nomad.solution$solution)
				degree <- K.opt[1:num.x]
				segments <- K.opt[(num.x+1):(2*num.x)]
				if(!is.null(z) ) I.opt <- K.opt[(2*num.x+1):(2*num.x+num.z)]
		}	
		else if(complexity=="degree") {
				degree <- as.integer(nomad.solution$solution[1:num.x])
				if(!is.null(z) ) {I.opt <- as.integer(nomad.solution$solution[(num.x+1):(num.x+num.z)])}
		}	
		else if(complexity=="knots")
		{
				segments <- as.integer(nomad.solution$solution[1:num.x])
				if(!is.null(z) ) {I.opt <- as.integer(nomad.solution$solution[(num.x+1):(num.x+num.z)])}
		}

		if(!is.null(z))
		{
				K.opt<-cbind(degree, segments , I.opt)
		}
		else {
				I.opt <- NULL
				K.opt<-cbind(degree, segments)
		}

		console <- printClear(console)
		console <- printPop(console)

    ## Set number of segments when degree==0 to 1 (or NA)

    segments[degree==0] <- 1

		if(any(degree==degree.max)) warning(paste(" optimal degree equals search maximum (", degree.max,"): rerun with larger degree.max",sep=""))
		if(any(segments==(segments.max+1))) warning(paste(" optimal segment equals search maximum (", segments.max+1,"): rerun with larger segments.max",sep=""))  

		basis.opt <- basis
		if(basis == "auto") basis.opt <- attributes(nomad.solution)$basis.opt
		
		cv.vec <- NULL
		basis.vec <- NULL
		KI.mat <- NULL

		crscv(K=K.opt,
					I=I.opt,
					basis=basis.opt,
					basis.vec=basis.vec,
					degree.max=degree.max, 
					segments.max=segments.max, 
					degree.min=degree.min, 
					segments.min=segments.min, 
					complexity=complexity,
					knots=knots,
					degree=degree,
					segments=segments,
					K.mat=KI.mat,
					restarts=NULL,
					lambda=NULL,
					lambda.mat=NULL,
					cv.objc=cv.min,
					cv.objc.vec=cv.vec,
					num.x=num.x,
					cv.func=cv.func)

}
