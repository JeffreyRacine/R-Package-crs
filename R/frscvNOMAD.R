## This function conducts factor regression spline
## cross-validation. It takes as input a data.frame xz containing a
## mix of numeric and factor predictors and a vector y. A range of
## arguments can be provided, and one can do search on both the degree
## and knots ("degree-knots") or the degree holding the number of
## knots (segments+1) constant or the number of knots (segments+1)
## holding the degree constant. Two basis types are supported
## ("additive", or "tensor") and the argument "auto" will choose the
## basis type automatically.

## Currently search is exhaustive taking basis.maxdim as the maximum
## number of the spline degree (0,1,...) and number of segments
## (1,2,...). This is a quadratic integer programming problem so
## ideally I require an IQP (MIQP for kernel-weighting)
## solver. Currently in R there is no such beast.

frscvNOMAD <- function(xz,
                       y,
                       basis.maxdim=5,
                       complexity=c("degree-knots","degree","knots"),
                       knots=c("quantiles","uniform"),
                       basis=c("additive","tensor","auto"),
                       cv.func=c("cv.ls","cv.gcv","cv.aic"),
                       degree=degree,
                       segments=segments, 
                       x0=x0, 
                       nmulti=0) {

		complexity <- match.arg(complexity)
		knots <- match.arg(knots)
		basis <- match.arg(basis)
		cv.func <- match.arg(cv.func)

		t1 <- Sys.time()

		cv.nomad <- function(x,
												 z=NULL,
												 y,
												 basis.maxdim, 
												 knots=knots,
												 basis=basis,
												 cv.func=cv.func, 
												 complexity=complexity, 
												 degree=degree, 
												 segments=segments, 
												 x0=x0, 
												 nmulti=nmulti) {

				if( missing(x) || missing(y) ||missing(basis.maxdim) ) stop(" you must provide input, x, y, and basis.maxdim")

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

				eval_cv <- function(input, params){
          
          complexity <- params$complexity
          segments <- params$segments
          degree <- params$degree
          x <- params$x
          y <- params$y
          z <- params$z
          knots <- params$knots
          cv.func <- params$cv.func
          basis <- params$basis
          
          num.x <- NCOL(x)
          
          if(complexity=="degree-knots") {
            K <- round(cbind(input[1:num.x],input[(num.x+1):(2*num.x)]))
          }	
          else if(complexity=="degree") {
            K<-round(cbind(input[1:num.x],segments))
          }	
          else if(complexity=="knots")
						{
              K<-round(cbind(degree, input[1:num.x]))
						}
          if(!is.null(z)) {
            num.z <- NCOL(z)
            I <- round(input[(2*num.x+1):(2*num.x+num.z)])
          } else {
            num.z <- 0
            I <- NULL
          }
          
          basis.opt <<-  basis;
          if(basis=="auto"){
            basis.opt<<-"additive"
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
              basis.opt <<-"tensor"
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
          
          return(cv)
				}
        
				##generate params
				params <- list()
				params$complexity <- complexity
				params$segments <- segments
				params$degree <- degree
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
				}

				if(complexity =="degree-knots") {
						if(is.null(x0))
								x0 <- c(sample(0:basis.maxdim,  num.x,  replace=T), sample(3:basis.maxdim, num.x, replace=T), sample(0:basis.maxdim, num.z, replace=T))
						bbin <-rep(1, num.x*2+num.z)
						#bounds segments cannot be samller than 2
						lb <- c(rep(0,num.x), rep(1, num.x), rep(0, num.z) )
						ub <- rep(basis.maxdim,num.x*2+num.z)

				}	
				else if(complexity=="degree") {
						if(is.null(x0))
								x0 <- c(sample(0:basis.maxdim,  num.x,  replace=T),  sample(0:basis.maxdim, num.z, replace=T))
						bbin <-rep(1, num.x+num.z)
						#bounds segments cannot be samller than 2
						lb <- c(rep(0,num.x),  rep(0, num.z) )
						ub <- rep(basis.maxdim,num.x+num.z)
				}	
				else if(complexity=="knots")
				{
						if(is.null(x0))
								x0 <- c(sample(1:basis.maxdim,  num.x,  replace=T),  sample(0:basis.maxdim, num.z, replace=T))
						bbin <-rep(1, num.x+num.z)
						#bounds segments cannot be samller than 2
						lb <- c(rep(1,num.x),  rep(0, num.z) )
						ub <- rep(basis.maxdim,num.x+num.z)
				}

				if(length(x0) != length(lb)) stop(" x0 and bounds have differing numbers of variables")

				#no constraints
				bbout <-c(0)

        ## This could be a passable parameter or set appropriately

				opts <-list("MAX_BB_EVAL"=500,
                    "MIN_MESH_SIZE"=1.0e-10,
                    "INITIAL_MESH_SIZE"=1.0e-02,
                    "MIN_POLL_SIZE"=1.0e-10)

#				opts <-list("MAX_BB_EVAL"=500,
#                    "MIN_MESH_SIZE"=0.00001,
#                    "INITIAL_MESH_SIZE"="0.1",
#                    "MIN_POLL_SIZE"=0.00001)

				solution<-snomadr(eval_f=eval_cv,
                          n=length(x0),
                          x0=as.numeric(x0),
                          bbin=bbin,
                          bbout=bbout,
                          lb=lb,
                          ub=ub,
                          nmulti=as.integer(nmulti),
                          opts=opts,
                          params=params);

		}

		console <- newLineConsole()
		console <- printPush("Solving by NOMAD...",console = console)

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
		}

		num.x <- ncol(x)
		n <- nrow(x)

		if(missing(x) || missing(y)) stop (" you must provide x and y")

		if(complexity=="degree") {
				if(missing(segments)) stop("segments missing for cross-validation of spline degree")
				if(length(segments)!=num.x) stop(" segments vector must be the same length as x")  
		} else if(complexity=="knots") {
				if(missing(degree)) stop("degree missing for cross-validation of number of spline knots")
				if(length(degree)!=num.x) stop(" degree vector must be the same length as x")
		}

		## For factor regression spline, if there is only one predictor
		## (i.e. num.x + num.z = 1) disable auto, set to additive (tensor in
		## this case, so don't waste time doing both).

		if(num.x+num.z==1 & basis == "auto") basis <- "additive"

		if(basis.maxdim < 1) stop(" basis.maxdim must be greater than or equal to 1")

		## we use a globle variable to store basis.opt
		## sovle by NOMAD

		nomad.solution <- cv.nomad(x,
															 z,
															 y,
															 basis.maxdim, 
															 knots=knots,
															 basis=basis,
															 cv.func=cv.func, 
															 complexity=complexity, 
															 degree=degree, 
															 segments=segments, 
															 x0=x0, 
															 nmulti=nmulti) 

		t2 <- Sys.time()

		##output

		cv.min <- nomad.solution$objective
		if(complexity=="degree-knots") {
				K.opt <- as.integer(nomad.solution$solution)
				degree <- K.opt[1:num.x]
				segments <- K.opt[(num.x+1):(2*num.x)]
				if(!is.null(z)) I.opt <- K.opt[(2*num.x+1):(2*num.x+num.z)]

		}	
		else if(complexity=="degree") {
				degree <- as.integer(nomad.solution$solution[1:num.x])
				if(!is.null(z)) {I.opt <- as.integer(nomad.solution$solution[(num.x+1):(num.x+num.z)])
				}

		}	
		else if(complexity=="knots")
		{
				segments <- as.integer(nomad.solution$solution[1:num.x])
				if(!is.null(z)) {I.opt <- as.integer(nomad.solution$solution[(num.x+1):(num.x+num.z)])
				}
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

		if(any(degree==basis.maxdim)) warning(paste(" optimal degree equals search maximum (", basis.maxdim,"): rerun with larger basis.maxdim",sep=""))
		if(any(segments==(basis.maxdim+1))) warning(paste(" optimal segment equals search maximum (", basis.maxdim+1,"): rerun with larger basis.maxdim",sep=""))  

		cv.vec <- NULL
		basis.vec <- NULL
		KI.mat <- NULL

		crscv(K=K.opt,
					I=I.opt,
					basis=basis.opt,
					basis.vec=basis.vec,
					basis.maxdim=basis.maxdim,
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
