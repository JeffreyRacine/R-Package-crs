# This function conducts kernel regression spline
## cross-validation. It takes as input a data.frame xz containing a
## mix of numeric and factor predictors and a vector y. A range of
## arguments can be provided, and one can do search on the bandwidths
## and both the degree and knots ("degree-knots") or the degree
## holding the number of knots (segments+1) constant or the number of
## knots (segments+1) holding the degree constant. Two basis types are
## supported ("additive" or "tensor") and the argument "auto" will
## choose the basis type automatically.

## Currently search is exhaustive taking basis.maxdim as the maximum
## number of the spline degree (0,1,...) and number of segments
## (1,2,...). This is a quadratic integer programming problem so
## ideally I require an IQP (MIQP for kernel-weighting)
## solver. Currently in R there is no such beast.

krscvNOMAD <- function(xz,
                       y,
                       basis.maxdim=5,
                       kernel.type=c("nominal","ordinal"),
                       restarts=0,
                       complexity=c("degree-knots","degree","knots"),
                       knots=c("quantiles","uniform"),
                       basis=c("additive","tensor","auto"),
                       cv.func=c("cv.ls","cv.gcv","cv.aic"),
                       degree=degree,
                       segments=segments, 
                       x0 = x0, 
                       nb_mads_runs=0) {

		complexity <- match.arg(complexity)
		knots <- match.arg(knots)
		basis <- match.arg(basis)
		cv.func <- match.arg(cv.func)  

		kernel.type <- match.arg(kernel.type)

		t1 <- Sys.time()

		cv.nomad <- function(x,
												 y,
												 z,
												 basis.maxdim,
												 restart,
												 num.restarts,
												 z.unique,
												 ind,
												 ind.vals,
												 nrow.z.unique,
												 kernel.type,
												 complexity=complexity,
												 knots=knots,
												 basis=basis,
												 cv.func=cv.func, 
												 x0=x0, 
												 nb_mads_runs=nb_mads_runs) {

				if( missing(x) || missing(y) ||missing(basis.maxdim) ) stop(" you must provide input, x, y, and basis.maxdim")

				## Presumes x (continuous predictors) exist, but z
				## (ordinal/nominal factors) can be optional

				n <- length(y)
				num.x <- NCOL(x)

				if(NROW(y) != NROW(x)) stop(" x and y have differing numbers of observations")

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
						z.unique <- params$z.unique
						ind <- params$ind
						ind.vals <- params$ind.vals
						nrow.z.unique <- params$nrow.z.unique
						kernel.type <- params$kernel.type

						num.x <- NCOL(x)
						num.z <- NCOL(z)

						if(complexity=="degree-knots") {
								K <- round(cbind(input[1:num.x],input[(num.x+1):(2*num.x)]))
								lambda <- input[(2*num.x+1):(2*num.x+num.z)]
						}	
						else if(complexity=="degree") {
								K<-round(cbind(input[1:num.x],segments))
								lambda <- input[(num.x+1):(num.x+num.z)]
						}	
						else if(complexity=="knots")
						{
								K<-round(cbind(degree, input[1:num.x]))
								lambda <- input[(num.x+1):(num.x+num.z)]
						}

						## When using weights= lambda of zero fails. Trivial to trap.
						## This will be a problem because we cannot change "input". I'll think it.
						## If cv=infinity, we will check in snomadr. 
						lambda <- ifelse(lambda <= 0, .Machine$double.eps, lambda)

						basis.opt <<-  basis;
						if(basis=="auto"){
								basis.opt<<-"additive"
								cv <- cv.kernel.spline(x=x,
																			 y=y,
																			 z=z,
																			 K=K,
																			 lambda=lambda,
																			 z.unique=z.unique,
																			 ind=ind,
																			 ind.vals=ind.vals,
																			 nrow.z.unique=nrow.z.unique,
																			 kernel.type=kernel.type,
																			 knots=knots,
																			 basis=basis.opt,
																			 cv.func=cv.func)
								cv.tensor <- cv.kernel.spline(x=x,
																							y=y,
																							z=z,
																							K=K,
																							lambda=lambda,
																							z.unique=z.unique,
																							ind=ind,
																							ind.vals=ind.vals,
																							nrow.z.unique=nrow.z.unique,
																							kernel.type=kernel.type,
																							knots=knots,
																							basis="tensor",
																							cv.func=cv.func)
								if(cv>cv.tensor){
										cv <- cv.tensor
										basis.opt <<-"tensor"
								}

						}
						else {
								cv <- cv.kernel.spline(x=x,
																			 y=y,
																			 z=z,
																			 K=K,
																			 lambda=lambda,
																			 z.unique=z.unique,
																			 ind=ind,
																			 ind.vals=ind.vals,
																			 nrow.z.unique=nrow.z.unique,
																			 kernel.type=kernel.type,
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
				params$z.unique <- z.unique
				params$ind <- ind
				params$ind.vals <- ind.vals
				params$nrow.z.unique <- nrow.z.unique
				params$kernel.type <- kernel.type
				# initial value
				num.z <- NCOL(z)

				if(complexity =="degree-knots") {
						if(is.null(x0))
								x0 <- c(sample(0:basis.maxdim,  num.x,  replace=T), sample(3:basis.maxdim, num.x, replace=T), runif(num.z))
						bbin <-c(rep(1, num.x*2),rep(0,  num.z))
						#bounds segments cannot be samller than 2
						lb <- c(rep(0,num.x), rep(1, num.x), rep(0, num.z) )
						ub <- c(rep(basis.maxdim,num.x*2), rep(1, num.z))

				}	
				else if(complexity=="degree") {
						if(is.null(x0))
								x0 <- c(sample(0:basis.maxdim,  num.x,  replace=T),  runif(num.z))
						bbin <-c(rep(1, num.x),rep(0,  num.z))
						lb <- c(rep(0,num.x),  rep(0, num.z) )
						ub <- c(rep(basis.maxdim,num.x), rep(1, num.z))
				}	
				else if(complexity=="knots")
				{
						if(is.null(x0))
								x0 <- c(sample(1:basis.maxdim,  num.x,  replace=T),  runif( num.z))
						bbin <-c(rep(1, num.x),rep(0,  num.z))
						#bounds segments cannot be samller than 2
						lb <- c(rep(1,num.x),  rep(0, num.z) )
						ub <- c(rep(basis.maxdim,num.x), rep(1, num.z))
				}

				if(length(x0) != length(lb)) stop(" x0 and bounds have differing numbers of variables")

				#no constraints
				bbout <-c(0)

				opts <-list("MAX_BB_EVAL"=500, "MIN_MESH_SIZE"=0.00001, "INITIAL_MESH_SIZE"="0.1", "MIN_POLL_SIZE"=0.00001)

				solution<-snomadr(eval_f=eval_cv, n=length(x0),x0=as.numeric(x0), bbin=bbin, bbout=bbout, lb=lb, ub=ub,nb_mads_runs=as.integer(nb_mads_runs),  opts=opts, params=params);

		}

		xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
		x <- xztmp$x
		z <- xztmp$z
		if(is.null(z)) stop(" categorical kernel smoothing requires ordinal/nominal predictors")

		z <- as.matrix(xztmp$z)
		num.z <- NCOL(z)
		z.unique <- uniquecombs(z)
		ind <-  attr(z.unique,"index")
		ind.vals <-  unique(ind)
		nrow.z.unique <- NROW(z.unique)
		num.x <- NCOL(x)
		n <- NROW(x)

		if(complexity=="degree") {
				if(missing(segments)) stop("segments missing for cross-validation of spline degree")
				if(length(segments)!=num.x) stop(" segments vector must be the same length as x")  
		} else if(complexity=="knots") {
				if(missing(degree)) stop("degree missing for cross-validation of number of spline knots")
				if(length(degree)!=num.x) stop(" degree vector must be the same length as x")
		}

		## For kernel regression spline, if there is only one continuous
		## predictor (i.e. num.x==1) disable auto, set to additive (which is
		## tensor in this case, so don't waste time doing both).

		if(num.x==1 & basis == "auto") basis <- "additive"

		if(basis.maxdim < 1) stop(" basis.maxdim must be greater than or equal to 1")

		console <- newLineConsole()
		console <- printPush("Solving by NOMAD...",console = console)

		basis.opt<<- "additive"

		## solve by NOMAD

		nomad.solution<-cv.nomad(x,
                             y,
                             z,
                             basis.maxdim=basis.maxdim,
                             restart=restart,
                             num.restarts=num.restarts,
                             z.unique=z.unique,
                             ind=ind,
                             ind.vals=ind.vals,
                             nrow.z.unique=nrow.z.unique,
                             kernel.type=kernel.type, 
                             complexity=complexity,
                             knots=knots,
                             basis=basis,
                             cv.func=cv.func, 
                             x0=x0, 
                             nb_mads_runs=nb_mads_runs) 
    
		t2 <- Sys.time()

		##output

		cv.min <- nomad.solution$objective
		if(complexity=="degree-knots") {
				K.opt <- as.integer(nomad.solution$solution[1:(2*num.x)])
				lambda.opt <- as.numeric(nomad.solution$solution[(2*num.x+1):(2*num.x+num.z)])
				degree <- K.opt[1:num.x]
				segments <- K.opt[(num.x+1):(2*num.x)]

		}	
		else if(complexity=="degree") {
				degree <- as.integer(nomad.solution$solution[1:num.x])
				lambda.opt <- as.numeric(nomad.solution$solution[(num.x+1):(num.x+num.z)])
				K.opt <-cbind(degree, segments)

		}	
		else if(complexity=="knots")
		{
				segments <- as.integer(nomad.solution$solution[1:num.x])
				lambda.opt <- as.numeric(nomad.solution$solution[(num.x+1):(num.x+num.z)])
				K.opt <-cbind(degree, segments)
		}

		console <- printClear(console)
		console <- printPop(console)

    print("here we are")
    print(nomad.solution$solution)
    print(segments)
    print(degree)
    print(basis.maxdim)

	if(any(degree==basis.maxdim)) warning(paste(" optimal degree equals search maximum (", basis.maxdim,"): rerun with larger basis.maxdim",sep=""))
	if(any(segments==(basis.maxdim+1))) warning(paste(" optimal segment equals search maximum (", basis.maxdim+1,"): rerun with larger basis.maxdim",sep=""))  

  ## We do not use the following parameters, so there may be errors if
  ## other functions will use them.
    
		cv.vec <- NULL
		lambda.mat <- NULL
		basis.vec <- NULL
		K.mat <- NULL  #We should check K.mat

		crscv(K=K.opt,
					I=NULL,
					basis=basis.opt,
					basis.vec=basis.vec,
					basis.maxdim=basis.maxdim,
					complexity=complexity,
					knots=knots,
					degree=degree,
					segments=segments,
					restarts=restarts,
					K.mat=K.mat,
					lambda=lambda.opt,
					lambda.mat=lambda.mat,
					cv.objc=cv.min,
					cv.objc.vec=cv.vec,
					num.x=num.x,
					cv.func=cv.func)

}
