crscv <- function(K,
                  I,
                  tensor,
                  tensor.vec,
                  max.K,
                  restarts,
                  K.mat,
                  lambda,
                  lambda.mat,
                  cv.func,
                  cv.func.vec) {
  
  tregcv = list(K=K,
                I=I,
                tensor=tensor,
                tensor.vec=tensor.vec,    
                max.K=max.K,
                restarts=restarts,
                K.mat=K.mat,
                lambda=lambda,
                lambda.mat=lambda.mat,
                cv.func=cv.func,
                cv.func.vec=cv.func.vec)

  class(tregcv) <- "crscv"

  tregcv
}

print.crscv <- function(x, ...){
  if(!is.null(x$lambda)&&is.null(x$I)) {
    cat("\nCategorical Regression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective Function Value : ",format(x$cv.func),sep=""),sep="")
    for(j in 1:length(x$K))
      cat(paste("\nOptimal spline degree for x[", j, "]: ",format(x$K[j]),sep=""),sep="")
    for(j in 1:length(x$lambda))
      cat(paste("\nOptimal bandwidth for z[", j, "]: ",format(x$lambda[j]),sep=""),sep="")
    cat(paste("\n\nMaximum Spline Degree for Search ", sQuote("max.K"), ": ",format(x$max.K),sep=""),sep="")
    cat(paste("\nTensor = ", x$tensor,sep=""))
    if(x$restarts>0) cat(paste("\nNumber of Restarts = ", format(x$restarts),sep=""),sep="")    
    cat("\n\n")
  } else if(!is.null(x$I)) {
    cat("\nFactor Regression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective Function Value : ",format(x$cv.func),sep=""),sep="")
    for(j in 1:length(x$K))
      cat(paste("\nOptimal spline degree for x[", j, "]: ",format(x$K[j]),sep=""),sep="")
    for(j in 1:length(x$I))
      cat(paste("\nInclusion for z[", j, "]: ",format(x$I[j]),sep=""),sep="")
    cat(paste("\n\nMaximum Spline Degree for Search ", sQuote("max.K"), ": ",format(x$max.K),sep=""),sep="")
    cat(paste("\nTensor = ", x$tensor,sep=""))
    cat("\n\n")
  } else {
    cat("\nRegression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective Function Value : ",format(x$cv.func),sep=""),sep="")
    for(j in 1:length(x$K))
      cat(paste("\nOptimal spline degree for x[", j, "]: ",format(x$K[j]),sep=""),sep="")
    cat(paste("\n\nMaximum Spline Degree for Search ", sQuote("max.K"), ": ",format(x$max.K),sep=""),sep="")
    cat(paste("\nTensor = ", x$tensor,sep=""))
    cat("\n\n")
  }
}

summary.crscv <- function(object, ...){
  print(object)
}
