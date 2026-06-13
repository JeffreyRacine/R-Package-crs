## This function conducts factor regression spline cross-validation
## using NOMAD. It takes as input a data.frame xz containing a mix of
## numeric and factor predictors and a vector y. A range of arguments
## can be provided, and one can do search on both the degree and knots
## ("degree-knots") or the degree holding the number of knots
## (segments+1) constant or the number of knots (segments+1) holding
## the degree constant. Three basis types are supported ("additive",
## "glp" or "tensor") and the argument "auto" will choose the basis
## type automatically.


frscvNOMAD <- function(xz,
                       y,
                       basis=c("additive","tensor","glp","auto"),
                       complexity=c("degree-knots","degree","knots"),
                       cv.df.min=1,
                       cv.func=c("cv.ls","cv.gcv","cv.aic"),
                       degree=degree,
                       degree.max=10,
                       degree.min=0,
                       display.nomad.progress=TRUE,
                       display.warnings=TRUE,
                       include=include,
                       initial.mesh.size.integer="1",
                       knots=c("quantiles","uniform", "auto"),
                       max.bb.eval=10000,
                       max.eval=NULL,
                       min.mesh.size.integer="1", #paste("r",sqrt(.Machine$double.eps),sep=""),
                       min.frame.size.integer="1",  #paste("r",sqrt(.Machine$double.eps),sep=""),
                       nmulti=0,
                       opts=list(),
                       random.seed=42,
                       segments=segments,
                       segments.max=10,
                       segments.min=1,
                       singular.ok=FALSE,
                       tau=NULL,
                       weights=NULL) {

  min.frame.size.integer.supplied <- !missing(min.frame.size.integer)

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)
  cv.func <- match.arg(cv.func)
  cv.maxPenalty <- resolve_cv_maxPenalty(NULL, y, weights = weights, cv.func = cv.func)

  if ( missing(include) || is.null(include)) {
    include <- NULL
  }
  if(degree.min < 0 ) degree.min <- 0
  if(segments.min < 1 ) segments.min <- 1
  if(degree.max < degree.min) degree.max <- (degree.min + 1)
  if(segments.max < segments.min) segments.max <- (segments.min + 1)
  if(missing(degree)) degree <- NULL
  if(missing(segments)) segments <- NULL

  ## Set DISPLAY_DEGREE to 0 if crs.messages=FALSE and
  ## DISPLAY_DEGREE is not provided

  if(!isTRUE(getOption("crs.messages")) && is.null(opts[["DISPLAY_DEGREE"]])) opts$"DISPLAY_DEGREE"=0
  opts <- merge.nomad4.fr.defaults(opts)

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
                       print.output=print.output,
                       nmulti=nmulti,
                       cv.df.min=cv.df.min,
                       tau=tau,
                       weights=weights,
                       singular.ok=singular.ok) {

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

    progress.env <- new.env(parent = emptyenv())
    progress.env$emit <- isTRUE(print.output)
    progress.env$progress.status <- .crs_progress_status_begin(
      enabled = isTRUE(print.output),
      surface = "nomad"
    )
    progress.env$nomad.progress <- NULL

    on.exit({
      .crs_progress_status_clear(progress.env$progress.status)
      progress.env$nomad.progress <- NULL
    }, add = TRUE)

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
      tau <- params$tau
      weights <- params$weights
      singular.ok <- params$singular.ok
      cv.df.min <- params$cv.df.min

      num.x <- NCOL(x)
      num.z <- 0
      I <- include

      if(complexity=="degree-knots") {
        K <- round(cbind(input[seq_len(num.x)], input[.crs_index_block(num.x, num.x)]))
        if(!is.null(z)  ) {
          num.z <- NCOL(z)
          I <- round(input[.crs_index_block(2 * num.x, num.z)])
        }
      }
      else if(complexity=="degree") {
        K<-round(cbind(input[seq_len(num.x)],segments))
        if(!is.null(z)  ) {
          num.z <- NCOL(z)
          I <- round(input[.crs_index_block(num.x, num.z)])
        }
      }
      else if(complexity=="knots")
      {
        K<-round(cbind(degree, input[seq_len(num.x)]))
        if(!is.null(z)  ) {
          num.z <- NCOL(z)
          I <- round(input[.crs_index_block(num.x, num.z)])
        }
      }

      basis.opt <-  basis;
      if(basis=="auto"){
        basis.opt <-"additive"

        cv <- cv.factor.spline.wrapper(x=x,
                                       y=y,
                                       z=z,
                                       K=K,
                                       I=I,
                                       knots=knots,
                                       basis=basis.opt,
                                       cv.func=cv.func,
                                       cv.df.min=cv.df.min,
                                       tau=tau,
                                       weights=weights,
                                       singular.ok=singular.ok)

        cv.tensor <- cv.factor.spline.wrapper(x=x,
                                              y=y,
                                              z=z,
                                              K=K,
                                              I=I,
                                              knots=knots,
                                              basis="tensor",
                                              cv.func=cv.func,
                                              cv.df.min=cv.df.min,
                                              tau=tau,
                                              weights=weights,
                                              singular.ok=singular.ok)
        if(cv > cv.tensor){
          cv <- cv.tensor
          basis.opt <-"tensor"
        }

        cv.glp <- cv.factor.spline.wrapper(x=x,
                                           y=y,
                                           z=z,
                                           K=K,
                                           I=I,
                                           knots=knots,
                                           basis="glp",
                                           cv.func=cv.func,
                                           cv.df.min=cv.df.min,
                                           tau=tau,
                                           weights=weights,
                                           singular.ok=singular.ok)
        if(cv>cv.glp){
          cv <- cv.glp
          basis.opt <-"glp"
        }

      } else {
        cv <- cv.factor.spline.wrapper(x=x,
                                       y=y,
                                       z=z,
                                       K=K,
                                       I=I,
                                       knots=knots,
                                       basis=basis.opt,
                                       cv.func=cv.func,
                                       cv.df.min=cv.df.min,
                                       tau=tau,
                                       weights=weights,
                                       singular.ok=singular.ok)
      }

      attr(cv, "basis.opt")<-basis.opt

      progress.env <- params$progress.env
      if (isTRUE(progress.env$emit)) {
        .crs_nomad_progress_maybe_restart(
          progress = progress.env$nomad.progress,
          input = input,
          degree = K[, 1L],
          segments = K[, 2L],
          aux = I
        )
        .crs_nomad_progress_eval(
          progress = progress.env$nomad.progress,
          degree = K[, 1L],
          segments = K[, 2L],
          aux = I,
          value = cv
        )
      }

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
    params$tau <- tau
    params$weights <- weights
    params$singular.ok <- singular.ok
    params$cv.df.min <- cv.df.min
    params$progress.env <- progress.env

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

    ## Save seed prior to setting

    seed.state <- .crs_capture_seed()

    set.seed(random.seed)

    if(is.null(xdegree)) xdegree <- rep(1,num.x) ## sample(degree.min:degree.max, num.x, replace=TRUE)
    if(is.null(xsegments)) xsegments <- rep(1,num.x) ## sample(segments.min:segments.max, num.x, replace=TRUE)
    if(is.null(xinclude)) xinclude <- sample(0:1, num.z, replace=TRUE)

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
      lb <- c(rep(degree.min, num.x), rep(0, num.z))
      ub <- c(rep(degree.max, num.x), rep(1, num.z))
    }
    else if(complexity=="knots")
    {
      x0 <- c(xsegments, xinclude)
      bbin <-rep(1, num.x+num.z)
      #bounds segments cannot be samller than 2
      lb <- c(rep(segments.min, num.x), rep(0, num.z))
      ub <- c(rep(segments.max, num.x), rep(1, num.z))
    }

    ## Test for ill-conditioned spline degree, reset upper bound
    ## if binding and start values.

    ill.conditioned <- check.max.spline.degree(x,rep(degree.max,num.x),display.warnings=display.warnings)
    degree.max.vec <- attr(ill.conditioned, "degree.max.vec")

    if(complexity != "knots") {
      ub[seq_len(num.x)] <- pmin(ub[seq_len(num.x)], degree.max.vec)
      x0[seq_len(num.x)] <- pmin(x0[seq_len(num.x)], degree.max.vec)
    }

    if(length(x0) != length(lb)) stop(" x0 and bounds have differing numbers of variables")

    #no constraints
    bbout <-c(0)

    ## Manual says precede by r means relative to up and lb... not
    ## quite what I was looking for

    x0.starts <- if (as.integer(nmulti) > 1L) {
      .crs_nomad_capture_start_matrix(
        x0 = x0,
        nstart = if (nmulti > 0L) as.integer(nmulti) else 1L,
        bbin = bbin,
        bbout = bbout,
        lb = lb,
        ub = ub,
        random.seed = random.seed,
        opts = opts
      )
    } else {
      NULL
    }

    progress.env$nomad.progress <- .crs_nomad_progress_begin(
      progress.status = progress.env$progress.status,
      label = "Selecting spline model",
      starts = if (!is.null(x0.starts)) x0.starts else matrix(as.numeric(x0), nrow = 1L, byrow = TRUE),
      aux.label = if (!is.null(z)) "include" else NULL,
      aux.type = "int"
    )

    solution <- if (!is.null(x0.starts)) {
      .crs_nomad_restart_sweep(eval.f=eval.cv,
                               n=length(x0),
                               starts=x0.starts,
                               bbin=bbin,
                               bbout=bbout,
                               lb=lb,
                               ub=ub,
                               random.seed=random.seed,
                               opts=opts,
                               display.nomad.progress=print.output,
                               params=params,
                               use.cache=FALSE)
    } else {
      snomadr(eval.f=eval.cv,
              n=length(x0),
              x0=as.numeric(x0),
              bbin=bbin,
              bbout=bbout,
              lb=lb,
              ub=ub,
              nmulti=as.integer(nmulti),
              random.seed=random.seed,
              opts=opts,
              display.nomad.progress=print.output,
              params=params)
    }

    .crs_nomad_progress_finish(progress.env$nomad.progress)
    progress.env$emit <- FALSE

    if(basis == "auto") {
      cv.basis <- eval.cv(solution$solution, params)
      attr(solution, "basis.opt") <- attributes(cv.basis)$basis.opt
      if(knots == "auto")
        attr(solution, "knots.opt") <- attributes(cv.basis)$knots.opt
    }
    else if(knots == "auto")
      attr(solution, "knots.opt") <- attributes(eval.cv(solution$solution, params))$knots.opt

    solution$degree.max.vec <- degree.max.vec

    ## Restore seed

    .crs_restore_seed(seed.state)

    return(solution)
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
      if(display.warnings) warning(paste(" the length of include (", length(include),") is not the same as the length of z (", num.z, ")",sep=""))
      include <- NULL
    }
    if (any(include <0) || any(include > 1)) include <- NULL
  }
  is.ordered.z <- xztmp$is.ordered.z
  num.x <- ncol(x)
  n <- nrow(x)

  if(missing(x) || missing(y)) stop (" you must provide x and y")

  if(complexity=="degree") {
    if(missing(segments)||is.null(segments)) stop(" segments missing for cross-validation of spline degree")
    if(length(segments)!=num.x) stop(" segments vector must be the same length as x")
    if(!is.null(degree) && length(degree) == num.x) { #check initial values should be in the bounds
      if(any(degree < degree.min)||any(degree>degree.max)) {
        if(display.warnings) warning(paste(" The provided initial values for the degree are not in the bounds.", sep=""))
        degree <- NULL
      }
    }
    else
      degree <- NULL
  } else if(complexity=="knots") {
    if(missing(degree)||is.null(degree)) stop(" degree missing for cross-validation of number of spline knots")
    if(length(degree)!=num.x) stop(" degree vector must be the same length as x")
    if(!is.null(segments) && length(segments) == num.x) {
      if(any(segments < segments.min)||any(segments > segments.max)) {
        if(display.warnings) warning(paste(" The provided initial values for the segments are not in the bounds.", sep=""))
        segments <- NULL
      }
    }
    else
      segments <- NULL
  }
  else {
    if(!is.null(degree) && length(degree) == num.x) { #check initial values should be in the bounds
      if(any(degree < degree.min)||any(degree>degree.max)) {
        if(display.warnings) warning(paste(" The provided initial values for the degree are not in the bounds.", sep=""))
        degree <- NULL
      }
    }
    else {
      degree <- NULL
    }
    if(!is.null(segments) && length(segments) == num.x) {
      if(any(segments < segments.min)||any(segments > segments.max)) {
        if(display.warnings) warning(paste(" The provided initial values for the segments are not in the bounds.", sep=""))
        segments <- NULL
      }
    }
    else
      segments <- NULL
  }

  geometry.roles <- if(complexity == "degree-knots") {
    rep("integer", 2 * num.x + num.z)
  } else {
    rep("integer", num.x + num.z)
  }

  opts$"EPSILON" <- .Machine$double.eps
  opts <- .crs_nomad_apply_eval_budget(opts,
                                       max.bb.eval = max.bb.eval,
                                       max.eval = max.eval,
                                       default.max.eval = 1000,
                                       context = "frscvNOMAD")
  opts <- .crs_nomad_apply_source_geometry(
    opts,
    roles = geometry.roles,
    initial.mesh.size.integer = initial.mesh.size.integer,
    min.mesh.size.integer = min.mesh.size.integer,
    min.frame.size.integer = if (min.frame.size.integer.supplied) min.frame.size.integer else NULL
  )

  if(display.nomad.progress) {
    if(!is.null(opts$DISPLAY_DEGREE)){
      if(opts$DISPLAY_DEGREE <= 0){
        display.nomad.progress <- FALSE
      }
    }
  }

  print.output <- display.nomad.progress

  ## For factor regression spline, if there is only one predictor
  ## (i.e. num.x + num.z = 1) disable auto, set to additive (tensor in
  ## this case, so don't waste time doing both).

  if((num.x+num.z==1) && (basis == "auto")) basis <- "additive"

  if(degree.max < 1 || segments.max < 1 ) stop(" degree.max or segments.max must be greater than or equal to 1")

  ## solve by NOMAD
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
                             print.output=display.nomad.progress,
                             nmulti=nmulti,
                             cv.df.min=cv.df.min,
                             tau=tau,
                             weights=weights,
                             singular.ok=singular.ok)

  t2 <- Sys.time()

  ##output

  I.opt <- NULL
  cv.min <- nomad.solution$objective
  if(isTRUE(all.equal(cv.min, cv.maxPenalty))) stop(" Search failed: restart with larger nmulti or smaller degree.max")
  basis.opt <- basis
  if(complexity=="degree-knots") {
    K.opt <- as.integer(nomad.solution$solution)
    degree <- K.opt[seq_len(num.x)]
    segments <- K.opt[.crs_index_block(num.x, num.x)]
    if(!is.null(z) ) I.opt <- K.opt[.crs_index_block(2 * num.x, num.z)]
  }
  else if(complexity=="degree") {
    degree <- as.integer(nomad.solution$solution[seq_len(num.x)])
    if(!is.null(z) ) {I.opt <- as.integer(nomad.solution$solution[.crs_index_block(num.x, num.z)])}
  }
  else if(complexity=="knots")
  {
    segments <- as.integer(nomad.solution$solution[seq_len(num.x)])
    if(!is.null(z) ) {I.opt <- as.integer(nomad.solution$solution[.crs_index_block(num.x, num.z)])}
  }

  if(!is.null(z))
  {
    K.opt<-cbind(degree, segments , I.opt)
  }
  else {
    I.opt <- NULL
    K.opt<-cbind(degree, segments)
  }

  ## Set number of segments when degree==0 to 1 (or NA)

  segments[degree==0] <- 1

  if(display.warnings) {
    if(any(degree==nomad.solution$degree.max.vec)) warning(paste(" optimal degree equals search maximum (", nomad.solution$degree.max.vec,"): rerun with larger degree.max",sep=""))
    if(any(segments==segments.max)) warning(paste(" optimal segment equals search maximum (", segments.max,"): rerun with larger segments.max",sep=""))
    if(!is.null(opts$MAX_BB_EVAL)){
      if(nmulti>0) {if(nmulti*opts$MAX_BB_EVAL <= nomad.solution$bbe) warning(paste(" MAX_BB_EVAL reached in NOMAD: perhaps use a larger value...", sep=""))}
      if(nmulti==0) {if(opts$MAX_BB_EVAL <= nomad.solution$bbe) warning(paste(" MAX_BB_EVAL reached in NOMAD: perhaps use a larger value...", sep="")) }
    }
  }

  basis.opt <- basis
  if(basis == "auto") basis.opt <- attributes(nomad.solution)$basis.opt

  knots.opt <- knots
  if(knots == "auto") knots.opt <- attributes(nomad.solution)$knots.opt

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
        knots=knots.opt,
        degree=degree,
        segments=segments,
        K.mat=KI.mat,
        restarts=nmulti,
        lambda=NULL,
        lambda.mat=NULL,
        cv.objc=cv.min,
        cv.objc.vec=cv.vec,
        num.x=num.x,
        cv.func=cv.func,
        tau=tau,
        nomad.restart.contract=nomad.solution$restart.contract,
        nomad.best.restart=nomad.solution$best.restart,
        nomad.restart.objectives=nomad.solution$restart.fval,
        nomad.restart.evaluations=nomad.solution$restart.results,
        nomad.summary=.crs_nomad_attach_effective_options(
          .crs_nomad_summary_from_solution(nomad.solution),
          opts
        ),
        cv.elapsed=as.numeric(difftime(t2, t1, units = "secs")))

}
