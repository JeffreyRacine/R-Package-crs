crscv <- function(K,
                  I,
                  basis,
                  basis.vec,
                  degree.max,
                  segments.max,
                  degree.min,
                  segments.min,
                  complexity,
                  knots,
                  degree,
                  segments,
                  restarts,
                  K.mat,
                  lambda,
                  lambda.mat,
                  cv.objc,
                  cv.objc.vec,
                  num.x,
                  cv.func,
                  tau,
                  nomad.restart.contract = NULL,
                  nomad.best.restart = NULL,
                  nomad.restart.objectives = NULL,
                  nomad.restart.evaluations = NULL,
                  nomad.summary = NULL) {

  tregcv = list(K=K,
                I=I,
                basis=basis,
                basis.vec=basis.vec,
                degree.max=degree.max,
                segments.max=segments.max,
                degree.min=degree.min,
                segments.min=segments.min,
                complexity=complexity,
                knots=knots,
                degree=degree,
                segments=segments,
                restarts=restarts,
                K.mat=K.mat,
                lambda=lambda,
                lambda.mat=lambda.mat,
                cv.objc=cv.objc,
                cv.objc.vec=cv.objc.vec,
                num.x=num.x,
                cv.func=cv.func,
                tau=tau)

  if (!is.null(nomad.restart.contract)) {
    tregcv$nomad.restart.contract <- nomad.restart.contract
  }
  if (!is.null(nomad.best.restart)) {
    tregcv$nomad.best.restart <- nomad.best.restart
  }
  if (!is.null(nomad.restart.objectives)) {
    tregcv$nomad.restart.objectives <- nomad.restart.objectives
  }
  if (!is.null(nomad.restart.evaluations)) {
    tregcv$nomad.restart.evaluations <- nomad.restart.evaluations
  }
  if (!is.null(nomad.summary)) {
    tregcv$nomad.summary <- nomad.summary
  }

  class(tregcv) <- "crscv"

  tregcv
}

print.crscv <- function(x, ...){

  if(!is.null(x$lambda)&&is.null(x$I)) {
    cat("\nCategorical Regression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective function: ", format(x$cv.func), sep=""))
    cat(paste("\nObjective function value: ",format(x$cv.objc),sep=""),sep="")

    cat(paste("\n\nKnot type: ", format(x$knots), sep=""))
    cat(paste("\nModel complexity proxy: ", format(x$complexity), sep=""))

    for(j in seq_along(x$degree))
      cat(paste("\nSpline degree/number of segments for x[", j, "]: ", format(x$degree[j]),"/",format(x$segments[j]),sep=""),sep="")
    if(!is.null(x$I)) for(j in seq_along(x$I))
      cat(paste("\nInclusion indicator for z[", j, "]: ",format(x$I[j]),sep=""),sep="")
    if(!is.null(x$lambda)) for(j in seq_along(x$lambda))
      cat(paste("\nBandwidth for  z[", j, "]: ",format(x$lambda[j]),sep=""),sep="")

    cat(paste("\n\nMaximum spline degree for search: ",format(x$degree.max),sep=""),sep="")
    cat(paste("\nBasis: ", x$basis,sep=""))
    if(x$restarts>0) cat(paste("\nNumber of restarts = ", format(x$restarts),sep=""),sep="")
    .crs_nomad_summary_print(x)
    cat("\n\n")
  } else if(!is.null(x$I)) {
    cat("\nFactor Regression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective function: ", format(x$cv.func), sep=""))
    cat(paste("\nObjective function value: ",format(x$cv.objc),sep=""),sep="")
    cat(paste("\n\nKnot type: ", format(x$knots), sep=""))
    cat(paste("\nModel complexity proxy: ", format(x$complexity), sep=""))

    for(j in seq_along(x$degree))
      cat(paste("\nSpline degree/number of segments for x[", j, "]: ", format(x$degree[j]),"/",format(x$segments[j]),sep=""),sep="")
    if(!is.null(x$I)) for(j in seq_along(x$I))
      cat(paste("\nInclusion indicator for z[", j, "]: ",format(x$I[j]),sep=""),sep="")
    if(!is.null(x$lambda)) for(j in seq_along(x$lambda))
      cat(paste("\nBandwidth for  z[", j, "]: ",format(x$lambda[j]),sep=""),sep="")

    cat(paste("\n\nMaximum spline degree for search: ",format(x$degree.max),sep=""),sep="")
    cat(paste("\nBasis: ", x$basis,sep=""))
    if(!is.null(x$restarts) && (x$restarts > 0)) cat(paste("\nNumber of restarts = ", format(x$restarts),sep=""),sep="")
    .crs_nomad_summary_print(x)
    cat("\n\n")
  } else {
    cat("\nRegression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective Function Value: ",format(x$cv.objc),sep=""),sep="")

    cat(paste("\n\nKnot type: ", format(x$knots), sep=""))
    cat(paste("\nModel complexity proxy: ", format(x$complexity), sep=""))

    for(j in seq_len(x$num.x))
      cat(paste("\nSpline degree/number of segments for x[", j, "]: ", format(x$degree[j]),"/",format(x$segments[j]),sep=""),sep="")

    if(!is.null(x$I)) for(j in seq_along(x$I))
      cat(paste("\nInclusion indicator for z[", j, "]: ",format(x$I[j]),sep=""),sep="")
    if(!is.null(x$lambda)) for(j in seq_along(x$lambda))
      cat(paste("\nBandwidth for  z[", j, "]: ",format(x$lambda[j]),sep=""),sep="")

    cat(paste("\n\nMaximum spline degree for search: ",format(x$degree.max),sep=""),sep="")
    cat(paste("\nBasis: ", x$basis,sep=""))
    if(!is.null(x$restarts) && (x$restarts > 0)) cat(paste("\nNumber of restarts = ", format(x$restarts),sep=""),sep="")
    .crs_nomad_summary_print(x)
    cat("\n\n")
  }
}

summary.crscv <- function(object, ...){
  print(object)
}
