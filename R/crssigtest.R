## (C) Jeffrey S. Racine (2011)

## This function performs an individual bootstrap test of significance
## for a crs regression model. If the index `i' is provided it
## computes the test for predictor `i', otherwise it computes the test
## for all predictors, one at a time.

## Two methods for imposing the null are implemented, a) a residual
## bootstrap and b) a reording of the predictor being tested (in
## place) that breaks any systematic relationship between the
## predictor and outcome.

crssigtest <- function(model = NULL,
                       index = NULL,
                       boot.num = 399,
                       boot.type = c("residual","reorder"),
                       random.seed = 42) {

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed <- TRUE
  } else {
    exists.seed <- FALSE
  }

  set.seed(random.seed)

  if(is.null(model)) stop(" you must provide a crs model")
  if(is.null(index)) index <- 1:NCOL(model$xz)
  if(index < 1 || index > NCOL(model$xz)) stop(" you must provide a valid index")
  P.vec <- numeric(length(index))
  F.vec <- numeric(length(index))  
  F.boot <- numeric(length=boot.num)
  F.boot.mat <- matrix(NA,nrow=boot.num,ncol=length(index))

  boot.type <- match.arg(boot.type)

  for(ii in 1:length(index)) {

    ## Get information from model each time the test is run

    model.degree <- model$degree
    model.lambda <- model$lambda
    model.segments <- model$segments
    model.basis <- model$basis
    model.knots <- model$knots

    ## Placeholders for null degree/lambda

    degree.restricted <- model.degree
    lambda.restricted <- model.lambda

    ## Determine whether variable(s) is numeric or factor

    xz.numeric <- FALSE
    if(is.numeric(model$xz[,index[ii]])) xz.numeric <- TRUE

    ## Test whether degree or bandwidth for the variable being tested
    ## is zero or one. First, the predictors can be in any order
    ## (numeric/factor) so we need to be careful to set the
    ## degree/lambda for the predictor being tested

    if(xz.numeric) {

      ## Predictor being tested is numeric, so need to determine
      ## degree index

      degree.index <- 0
      for(jj in 1:index[ii]) if(is.numeric(model$xz[,jj])) degree.index <- degree.index + 1

      ## If the degree is zero (i.e. cross-validation has determined a
      ## variable is `irrelevant'), allow the model to be included
      ## using a degree 1 fit

      model.degree[degree.index] <- ifelse(model$degree[degree.index]==0,1,model$degree[degree.index])

      ## Set the null degree to zero

      degree.restricted[degree.index] <- 0

    } else {

      ## Predictor being tested is categorical, so need to determine
      ## lambda index

      lambda.index <- 0
      for(jj in 1:index[ii]) if(!is.numeric(model$xz[,jj])) lambda.index <- lambda.index + 1

      ## If the bandwidth is one (i.e. cross-validation has determined
      ## a variable is `irrelevant'), allow the model to be included
      ## using the `frequency' fit (bandwidth of basically zero)

      model.lambda[lambda.index] <- ifelse(isTRUE(all.equal(model$lambda[lambda.index],1)),.Machine$double.eps,model$lambda[lambda.index])

      ## Set the null bandwidth to one

      lambda.restricted[lambda.index] <- 1

    }

    model.unrestricted <- crs(xz=model$xz,
                              y=model$y,
                              cv="none",
                              degree=model.degree,
                              lambda=model.lambda,
                              segments=model.segments,
                              basis=model.basis,
                              knots=model.knots)

    model.restricted <- crs(xz=model$xz,
                            y=model$y,
                            cv="none",
                            degree=degree.restricted,
                            lambda=lambda.restricted,
                            segments=model.segments,
                            basis=model.basis,
                            knots=model.knots)

    ## Compute the pseudo F-value under the null

    uss <- sum(residuals(model.unrestricted)^2)
    rss <- sum(residuals(model.restricted)^2)

    df1 <- model$nobs-model.unrestricted$k
    df2 <- model.unrestricted$k-model.restricted$k

    F.df <- df1/df2
  
    F.pseudo <- F.df*(rss-uss)/uss
    F.vec[ii] <- F.pseudo

    if(boot.type=="reorder") xz.boot <- model$xz

    for(b in 1:boot.num) {
    
      ## Bootstrap sample under the null and recompute the
      ## `restricted' and `unrestricted' models for the null sample

      if(boot.type=="reorder") {
      
        xz.boot[,index[ii]] <- sample(model$xz[,index[ii]],replace=T)

        model.unrestricted.boot <- crs(xz=xz.boot,
                                       y=model$y,
                                       cv="none",
                                       degree=model.degree,
                                       segments=model.segments,
                                       lambda=model.lambda,
                                       basis=model.basis,
                                       knots=model.knots)

        model.restricted.boot <- crs(xz=xz.boot,
                                     y=model$y,
                                     cv="none",
                                     degree=degree.restricted,
                                     segments=model.segments,
                                     lambda=lambda.restricted,
                                     basis=model.basis,
                                     knots=model.knots)

      } else {

        y.boot <- fitted(model.restricted) + sample(residuals(model.unrestricted),replace=T)

        model.unrestricted.boot <- crs(xz=model$xz,
                                       y=y.boot,
                                       cv="none",
                                       degree=model.degree,
                                       lambda=model.lambda,
                                       segments=model.segments,
                                       basis=model.basis,
                                       knots=model.knots)

        model.restricted.boot <- crs(xz=model$xz,
                                     y=y.boot,
                                     cv="none",
                                     degree=degree.restricted,
                                     lambda=lambda.restricted,
                                     segments=model.segments,
                                     basis=model.basis,
                                     knots=model.knots)

      }


      ## Recompute the pseudo F-value under the null
      
      uss.boot <- sum(residuals(model.unrestricted.boot)^2)
      
      rss.boot <- sum(residuals(model.restricted.boot)^2)    
      
      F.boot[b] <- F.df*(rss.boot-uss.boot)/uss.boot

    }

    F.boot.mat[,ii] <- F.boot

    P.vec[ii] <- mean(ifelse(F.boot > F.pseudo, 1, 0))

  }

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
  
  return(list(index=index,
              P=P.vec,
              F=F.vec,
              F.boot=F.boot.mat,
              df1=df1,
              df2=df2))

} 
