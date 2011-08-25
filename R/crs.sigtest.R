## This function performs an individual bootstrap test of significance
## for a crs regression model. If the index `i' is provided it
## computes the test for predictor `i', otherwise it computes the test
## for all predictors, one at a time.

crs.sigtest <- function(model=NULL,index=NULL,B=399) {

  if(is.null(model)) stop(" you must provide a crs model")
  if(is.null(index)) index <- 1:NCOL(model$xz)
  if(index < 1 || index > NCOL(model$xz)) stop(" you must provide a valid index")
  P.vec <- numeric(length(index))
  F.vec <- numeric(length(index))  
  F.boot <- numeric(length=B)  

  for(i in 1:length(index)) {

    ## Get information from model each time the test is run

    model.degree <- model$degree
    model.lambda <- model$lambda
    model.segments <- model$segments
    model.basis <- model$basis
    model.knots <- model$knots

    ## Placeholders for null degree/lambda

    degree.null <- model.degree
    lambda.null <- model.lambda

    ## Determine whether variable(s) is numeric or factor

    xz.numeric <- FALSE
    if(is.numeric(model$xz[,index[i]])) xz.numeric <- TRUE

    ## Test whether degree or bandwidth for the variable being tested
    ## is zero or one. First, the predictors can be in any order
    ## (numeric/factor) so we need to be careful to set the
    ## degree/lambda for the predictor being tested

    if(xz.numeric) {
      ## Predictor being tested is numeric, so need to determine
      ## degree index
      degree.index <- 0
      for(j in 1:index[i]) if(is.numeric(model$xz[,j])) degree.index <- degree.index + 1
      ## If the degree is zero (i.e. cross-validation has determined a
      ## variable is `irrelevant'), allow the model to be included
      ## using a degree 1 fit
      model.degree[degree.index] <- ifelse(model$degree[degree.index]==0,1,model$degree[degree.index])
      ## Set the null degree to zero
      degree.null[degree.index] <- 0
    } else {
      ## Predictor being tested is categorical, so need to determine
      ## lambda index
      lambda.index <- 0
      for(j in 1:index[i]) if(!is.numeric(model$xz[,j])) lambda.index <- lambda.index + 1
      ## If the bandwidth is one (i.e. cross-validation has determined
      ## a variable is `irrelevant'), allow the model to be included
      ## using the `frequency' fit (bandwidth of basically zero)
      model.lambda[lambda.index] <- ifelse(isTRUE(all.equal(model$lambda[lambda.index],1)),.Machine$double.eps,model$lambda[lambda.index])
      ## Set the null bandwidth to one
      lambda.null[lambda.index] <- 1
    }

    model.unrestricted <- crs(xz=model$xz,y=model$y,cv="none",degree=model.degree,segments=model.segments,lambda=model.lambda,basis=model.basis,knots=model.knots)
    model.null <- crs(xz=model$xz,y=model$y,cv="none",degree=degree.null,segments=model.segments,lambda=lambda.null,basis=model.basis,knots=model.knots)

    ## Compute the pseudo F-value under the null

    uss <- sum(residuals(model.unrestricted)^2)
    rss <- sum(residuals(model.null)^2)
  
    F.pseudo <- (rss-uss)/uss
    F.vec[i] <- F.pseudo

    for(b in 1:B) {
    
      ## Bootstrap sample under the null
      
      y.boot <- fitted(model.null) + sample(residuals(model.null),replace=T)
      
      ## Recompute the `restricted' and `unrestricted' models for the
      ## null sample
      
      model.unrestricted.boot <- crs(xz=model$xz,y=y.boot,cv="none",degree=model.degree,segments=model.segments,lambda=model.lambda,basis=model.basis,knots=model.knots)
      model.null.boot <- crs(xz=model$xz,y=y.boot,cv="none",degree=degree.null,segments=model.segments,lambda=lambda.null,basis=model.basis,knots=model.knots)
      
      ## Recompute the pseudo F-value under the null
      
      rss.boot <- sum(residuals(model.null.boot)^2)    
      uss.boot <- sum(residuals(model.unrestricted.boot)^2)
      
      F.boot[b] <- (rss.boot-uss.boot)/uss.boot
    
    }

    P.vec[i] <- mean(ifelse(F.boot > F.pseudo, 1, 0))

  }

  return(list(index=index,P=P.vec,F=F.vec))

} 
