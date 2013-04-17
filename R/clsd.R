## These functions are for (currently univariate) logspline density
## estimation written by racinej@mcmaster.ca (Jeffrey S. Racine). They
## make use of spline routines in the crs package (available on
## CRAN). The approach involves joint selection of the degree and
## knots in contrast to the typical approach (e.g. Kooperberg and
## Stone) that sets the degree to 3 and optimizes knots only. Though
## more computationally demanding, the estimators are more efficient
## on average.

integrate.trapezoidal <- function(x,y) {

  ## This function will compute the cumulative integral at each sample
  ## realization using the Newton-Cotes trapezoidal rule and the
  ## cumsum function as we need to compute this in a computationally
  ## efficient manner. It can be used to return the distribution
  ## function from the density function etc.

  n <- length(x)
  rank.x <- rank(x)
  order.x <- order(x)
  y <- y[order.x]
  x <- x[order.x]
  int.vec <- numeric(length(x))
  int.vec[2:n] <- cumsum((x[2:n] - x[2:n-1]) * (y[2:n] + y[2:n-1]) / 2)

  return(int.vec[rank.x])

}

par.init <- function(degree,segments) {

  ## This function initializes parameters for search along with upper
  ## and lower bounds if appropriate.

  dim.p <- degree+segments

  ## The weights for the linear tails must be non-positive. The lower
  ## bound places a maximum bound on how quickly the tails are allowed
  ## to die off. Trial and error suggests the values below seem to be
  ## appropriate for a wide range of (univariate)
  ## distributions. Kooperberg suggests that in order to get the
  ## constraint theta < 0 use theta <= -epsilon for some small epsilon
  ## > 0. We therefore use sqrt machine epsilon.

  par.ub <- - sqrt(.Machine$double.eps)
  par.init <- c(runif(1,-10,par.ub),rnorm(dim.p-2),runif(1,-10,par.ub))
  par.upper <- c(par.ub,rep(Inf,dim.p-2),par.ub)

  return(list(par.init=par.init,
              par.upper=par.upper))

}

clsd <- function(x=NULL,
                 beta=NULL,
                 xeval=NULL,
                 degree=NULL,
                 segments=NULL,
                 degree.min=2,
                 degree.max=5,
                 segments.min=1,
                 segments.max=5,
                 lbound=NULL,
                 ubound=NULL,
                 basis="tensor",
                 knots="quantiles",
                 penalty=c("aic","sic","cv","none"),
                 deriv.index=1,
                 deriv=0,
                 do.break=FALSE,
                 do.gradient=TRUE,
                 er=NULL,
                 monotone=TRUE,
                 n.integrate=1.0e+03,
                 nmulti=1,
                 method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
                 verbose=FALSE) {
  
  ptm <- system.time("")

  if(is.null(x)) stop(" You must provide data")

  ## If no er is provided use the following ad-hoc rule which attempts
  ## to ensure we cover the support of the variable for distributions
  ## with moments. This gets the chi-square, t, and Gaussian for n >=
  ## 100 with all degrees of freedom and df=1 is perhaps the worst
  ## case scenario. This rule delivers er = 0.43429448, 0.21714724,
  ## 0.14476483, 0.10857362, 0.08685890, and 0.0723824110, for n = 10,
  ## 10^2, 10^3, 10^4, 10^5, and 10^6. It is probably too aggressive
  ## for the larger samples but one can override - the code traps for
  ## non-finite integration and issues a message when this occurs
  ## along with a suggestion.

  if(!is.null(er) && er < 0) stop(" er must be non-negative")
  if(is.null(er)) er <- 1/log(length(x))

  penalty <- match.arg(penalty)
  method <- match.arg(method)

  fv <- NULL

  if(is.null(beta)) {

    ## If no parameters are provided presume intention is to run
    ## maximum likelihood estimation to obtain the parameter
    ## estimates.
    
    ptm <- ptm + system.time(ls.ml.out <- ls.ml(x=x,
                                                degree.min=degree.min,
                                                degree.max=degree.max,
                                                segments.min=segments.min,
                                                segments.max=segments.max,
                                                lbound=lbound,
                                                ubound=ubound,
                                                method=method,
                                                nmulti=nmulti,
                                                do.gradient=do.gradient,
                                                er=er,
                                                do.break=do.break,
                                                penalty=penalty,
                                                monotone=monotone,
                                                verbose=verbose))
    
    beta <- ls.ml.out$beta
    degree <- ls.ml.out$degree
    segments <- ls.ml.out$segments
    fv <- ls.ml.out$fv

  }

  if(is.null(degree)) stop(" You must provide spline degree")
  if(is.null(segments)) stop(" You must provide number of segments")

  ## To obtain the constant of integration for B-spline bases, we need
  ## to compute log(integral exp(P%*%beta)) so we take an equally
  ## spaced extended range grid of length n plus the sample
  ## realizations (min and max of sample therefore present for what
  ## follows), and evaluation points xeval if they exist.

  if(is.null(xeval)) {
    ## x will be the first 1:length(x) elements in object[rank.xnorm]
    er <- extendrange(x,f=er)
    if(!is.null(lbound)) er[1] <- lbound
    if(!is.null(ubound)) er[2] <- ubound    
    xnorm <- c(x,seq(er[1],er[2],length=n.integrate))
    rank.xnorm <- rank(xnorm)
    order.xnorm <- order(xnorm)
    xnorm <- xnorm[order.xnorm]
  } else {
    ## xeval will be the first 1:length(xeval) elements in
    ## object[rank.xnorm]
    er <- extendrange(c(x,xeval),f=er)
    if(!is.null(lbound)) er[1] <- lbound
    if(!is.null(ubound)) er[2] <- ubound    
    xnorm <- c(xeval,x,seq(er[1],er[2],length=n.integrate))
    rank.xnorm <- rank(xnorm)
    order.xnorm <- order(xnorm)
    xnorm <- xnorm[order.xnorm]

    if(min(xeval) < er[1] | max(xeval) > er[2]) warning(" evaluation data extends beyond the range of `er'")
  }

  ## Charles Kooperberg has a manuscript "Statistical Modeling with
  ## Spline Functions', Jan 5 2006 on his web page. Chapter 6, page
  ## 286, figure 6.7 reveals a hybrid spline basis that in essence
  ## appears to drop two columns from my B-spline with 2 segments
  ## added artificially. This has the effect of removing the two bases
  ## that were delivering weight in tails leading to `kinks'. Hat-tip
  ## to Charles for his clear descriptions. Note we require the same
  ## for the derivatives below. Note that this logspline basis does
  ## not have the B-spline property that the pointwise sum of the
  ## bases is 1 everywhere.

  ptm <- ptm + system.time(suppressWarnings(Pnorm <- prod.spline(x=x,
                                                                 xeval=xnorm,
                                                                 K=cbind(degree,segments+if(monotone){2}else{0}),
                                                                 knots=knots,
                                                                 basis=basis)))
  
  if(monotone) Pnorm <- Pnorm[,-c(2,degree+segments+1)]

  ## Compute the normalizing constant so that the estimate
  ## integrates to one. We append linear splines to the B-spline
  ## basis to generate exponentially declining tails (K=cbind(1,1)
  ## creates the linear basis).

  ptm <- ptm + system.time(suppressWarnings(P.lin <- prod.spline(x=x,
                                                                 xeval=xnorm,
                                                                 K=cbind(1,1),
                                                                 knots=knots,
                                                                 basis=basis)))
  
  ## We append the linear basis to the left and rightmost polynomial
  ## bases. We match the slope of the linear basis to that of the
  ## polynomial basis at xmin/xmax (note that
  ## Pnorm[xnorm==max(x),-ncol(Pnorm)] <- 0 is there because the
  ## gsl.bspline values at the right endpoint are very small but not
  ## exactly zero but want to rule out any potential issues hence
  ## set them correctly to zero)
  
  Pnorm[xnorm<min(x),] <- 0
  Pnorm[xnorm>max(x),] <- 0
  Pnorm[xnorm==max(x),-ncol(Pnorm)] <- 0
  P.left <- as.matrix(P.lin[,1])
  P.right <- as.matrix(P.lin[,2])
  
  ## We want the linear segment to have the same slope as the
  ## polynomial segment it connects with and to match at the joint
  ## hence conduct some carpentry at the left boundary.
  
  index <- which(xnorm==min(x))
  index.l <- index+1
  index.u <- index+5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.poly.left <- as.numeric((Pnorm[index.u,1]-Pnorm[index.l,1])/(x.u-x.l))
  index.l <- index+1
  index.u <- index+5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.linear.left <- as.numeric((P.left[index.u]-P.left[index.l])/(x.u-x.l))
  
  ## Complete carpentry at the right boundary.
  
  index <- which(xnorm==max(x))
  index.l <- index-1
  index.u <- index-5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.poly.right <- as.numeric((Pnorm[index.u,ncol(Pnorm)]-Pnorm[index.l,ncol(Pnorm)])/(x.u-x.l))
  index.l <- index-1
  index.u <- index-5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.linear.right <- as.numeric((P.right[index.u]-P.right[index.l])/(x.u-x.l))
  
  ## Here are the linear segments with matching slopes XXX patch up
  ## deriv outside range of data needed as well XXX
  
  P.left <- as.matrix(P.left-1)*slope.poly.left/slope.linear.left+1
  P.right <- as.matrix(P.right-1)*slope.poly.right/slope.linear.right+1
  
  P.left[xnorm>=min(x),1] <- 0
  P.right[xnorm<=max(x),1] <- 0
  
  Pnorm[,1] <- Pnorm[,1]+P.left
  Pnorm[,ncol(Pnorm)] <- Pnorm[,ncol(Pnorm)]+P.right
  
  if(ncol(Pnorm)!=length(beta)) stop(paste(" Incompatible arguments: beta must be of dimension ",ncol(Pnorm),sep=""))

  Pnorm.beta <- as.numeric(Pnorm%*%as.matrix(beta))

  if(deriv > 0) {

    ptm <- ptm + system.time(suppressWarnings(Pnorm.deriv <- prod.spline(x=x,
                                                                         xeval=xnorm,
                                                                         K=cbind(degree,segments+if(monotone){2}else{0}),
                                                                         knots=knots,
                                                                         basis=basis,
                                                                         deriv.index=deriv.index,
                                                                         deriv=deriv)))
    
    if(monotone) Pnorm.deriv <- Pnorm.deriv[,-c(2,degree+segments+1)]

    ptm <- ptm + system.time(suppressWarnings(P.lin <- prod.spline(x=x,
                                                                   xeval=xnorm,
                                                                   K=cbind(1,1),
                                                                   knots=knots,
                                                                   basis=basis,
                                                                   deriv.index=deriv.index,
                                                                   deriv=deriv)))
    
    ## For the derivative bases on the extended range `xnorm', above
    ## and below max(x)/min(x) we assign the bases to constants
    ## (zero). We append the linear basis to the left and right of the
    ## bases. The left basis takes on linear values to the left of
    ## min(x), zero elsewhere, the right zero to the left of max(x),
    ## linear elsewhere.
    
    Pnorm.deriv[xnorm<min(x),] <- 0
    Pnorm.deriv[xnorm>max(x),] <- 0
    P.left <- as.matrix(P.lin[,1])
    P.left[xnorm>=min(x),1] <- 0
    P.right <- as.matrix(P.lin[,2])
    P.right[xnorm<=max(x),1] <- 0
    Pnorm.deriv[,1] <- Pnorm.deriv[,1]+P.left
    Pnorm.deriv[,ncol(Pnorm.deriv)] <- Pnorm.deriv[,ncol(Pnorm.deriv)]+P.right
    P.deriv.beta <- as.numeric(Pnorm.deriv%*%beta)

  } else {

    f.deriv <- NULL

  }

  ## Compute the constant of integration to normalize the density
  ## estimate so that it integrates to one.

  norm.cumsum.integrate <- integrate.trapezoidal(xnorm,exp(Pnorm.beta))
  norm.constant <- norm.cumsum.integrate[length(xnorm)]
  log.norm.constant <- log(norm.constant)

  if(!is.finite(log.norm.constant))
    stop(" integration not finite - perhaps reduce lb for endpoint weights")

  ## For the distribution, compute the density over the extended
  ## range, then return values corresponding to either the sample x or
  ## evaluation x (xeval) based on integration over the extended range
  ## for the xnorm points (xnorm contains x and xeval - this ought to
  ## ensure integration to one).

  ## f.norm is the density evaluated on the extended range (including
  ## sample observations and evaluation points if the latter exist),
  ## F.norm the distribution evaluated on the extended range.

  f.norm <- exp(Pnorm.beta-log.norm.constant)
  F.norm <- integrate.trapezoidal(xnorm,f.norm)
  if(deriv>0) f.norm.deriv <- as.numeric(f.norm*P.deriv.beta)

  ## Next, strip off the values of the distribution corresponding to
  ## either sample x or evaluation xeval

  if(is.null(xeval)) {
    f <-   f.norm[rank.xnorm][1:length(x)]
    F <-   F.norm[rank.xnorm][1:length(x)]
    if(deriv>0) f.deriv <- f.norm.deriv[rank.xnorm][1:length(x)]
    P <-   Pnorm[rank.xnorm,][1:length(x),]
    P.beta <- Pnorm.beta[rank.xnorm][1:length(x)]
  } else {
    f <-   f.norm[rank.xnorm][1:length(xeval)]
    F <-   F.norm[rank.xnorm][1:length(xeval)]
    if(deriv>0) f.deriv <- f.norm.deriv[rank.xnorm][1:length(xeval)]
    P <-   Pnorm[rank.xnorm,][1:length(xeval),]
    P.beta <- Pnorm.beta[rank.xnorm][1:length(xeval)]
  }

  clsd.return <- list(density=f,
                      density.deriv=f.deriv,
                      distribution=F,
                      density.er=f.norm,
                      distribution.er=F.norm,
                      xer=xnorm,
                      Basis.beta=P.beta,
                      Basis.beta.er=Pnorm.beta,
                      P=P,
                      Per=Pnorm,
                      logl=sum(P.beta-log.norm.constant),## issue XXXX this is potentially for evaluation data
                      constant=norm.constant,
                      degree=degree,
                      segments=segments,
                      knots=knots,
                      basis=basis,
                      nobs=length(x),
                      beta=beta,
                      fv=fv,
                      er=er,
                      penalty=penalty,
                      nmulti=nmulti,
                      x=x,
                      ptm=ptm)
  
  class(clsd.return) <- "clsd"
  return(clsd.return)

}

sum.log.density <- function(beta,
                            x,
                            degree,
                            segments,
                            lbound=NULL,
                            ubound=NULL,
                            basis="tensor",
                            knots="quantiles",
                            er=NULL,
                            n.integrate=1.0e+03,
                            penalty=c("aic","sic","cv","none"),
                            monotone=TRUE) {

  penalty <- match.arg(penalty)
  if(missing(x)) stop(" You must provide data")
  if(missing(beta)) stop(" You must provide coefficients")
  if(missing(degree)) stop(" You must provide spline degree")
  if(missing(segments)) stop(" You must provide number of segments")

  output <- clsd(beta=beta,
                 x=x,
                 degree=degree,
                 segments=segments,
                 lbound=lbound,
                 ubound=ubound,
                 basis=basis,
                 knots=knots,
                 er=er,
                 n.integrate=n.integrate,
                 monotone=TRUE)

  logl <- output$logl
  f.hat <- output$density

  complexity <- degree+segments-3

  if(penalty=="aic") {

    ## AIC (note the penalty log(n) is Schwarz-Bayes which heavily
    ## penalizes overfitting, while the penalty 3 is used by
    ## Kooperberg & Stone (1991), while 2 corresponds to LSCV)

    ## In his R package Kooperberg (which cites his 1997 Annals paper)
    ## uses -2 * loglikelihood + penalty * (number of knots - 1) with
    ## default log(samplesize) as in BIC. There is no penalty for the
    ## cubic spline itself. My segments is his number of knots minus
    ## one.

    return(2*logl-2*complexity)

  } else if(penalty=="sic") {

    ## Schwarz-Bayes IC (The BIC is an asymptotic result derived under
    ## the assumptions that the data distribution is in an exponential
    ## family).  Note the penalty log(n) is Schwarz-Bayes. The penalty
    ## 3 is used by Kooperberg & Stone, while 2 corresponds to LSCV.

    return(2*logl-log(length(f.hat))*complexity)

  } else if(penalty=="cv") {

    ## For delete-one-cross-validation, to compute the diagonal of the
    ## hat matrix we invoke a ghost call to lm to get the LU
    ## decomposition etc. which simply computes
    ## P%*%solve(t(P)%*%P)%*%t(P)) to provide the effective number of
    ## parameters. Delete-one ML (call to lm.influence verified to
    ## produce diag[P%*%solve(t(P)%*%P)%*%t(P)]...)

    h <- lm.influence(lm(rep(1,nrow(output$P))~output$P))$hat
    k <- round(sum(h))

    ## Compute 1-h, check for case where h=1

    one.minus.h <- (1-ifelse(h < 1, h, h-.Machine$double.eps))

    ## Use the delete-one estimator for the ml-cv function

    return(sum(log(f.hat/one.minus.h)))

  } else {

    ## No penalty

    return(logl)

  }

}

sum.log.density.gradient <- function(beta,
                                     x,
                                     degree,
                                     segments,
                                     lbound=NULL,
                                     ubound=NULL,
                                     basis="tensor",
                                     knots="quantiles",
                                     er=NULL,
                                     n.integrate=1.0e+03,
                                     penalty=c("aic","sic","cv","none"),
                                     monotone=TRUE) {

  ## Note - gradients are only ever computed for the sample
  ## realizations used to compute the log likelihood function

  penalty <- match.arg(penalty)
  if(missing(x)) stop(" You must provide data")
  if(missing(beta)) stop(" You must provide coefficients")
  if(missing(degree)) stop(" You must provide spline degree")
  if(missing(segments)) stop(" You must provide number of segments")

  output <- clsd(beta=beta,
                 x=x,
                 degree=degree,
                 segments=segments,
                 lbound=lbound,
                 ubound=ubound,
                 basis=basis,
                 knots=knots,
                 er=er,
                 n.integrate=n.integrate,
                 monotone=TRUE)

  exp.P.beta.P <- exp(output$Basis.beta.er)*output$Per
  int.exp.P.beta.P <- numeric(length=ncol(output$Per))
  for(i in 1:ncol(output$P)) int.exp.P.beta.P[i] <- integrate.trapezoidal(output$xer,exp.P.beta.P[,i])[length(output$xer)]

  if(penalty=="aic"|penalty=="sic") {
    return(2*(colSums(output$P)-length(x)*int.exp.P.beta.P/output$constant))
  } else {
    return(colSums(output$P)-length(x)*int.exp.P.beta.P/output$constant)
  }

}

ls.ml <- function(x,
                  degree.min=2,
                  segments.min=1,
                  degree.max=5,
                  segments.max=5,
                  lbound=NULL,
                  ubound=NULL,
                  do.break=FALSE,
                  do.gradient=TRUE,
                  maxit=10^5,
                  nmulti=1,
                  er=NULL,
                  method=NULL,
                  n.integrate=1.0e+03,
                  basis="tensor",
                  knots="quantiles",
                  penalty=c("aic","sic","cv","none"),
                  monotone=TRUE,
                  verbose=FALSE,
                  max.attempts=25,
                  random.seed=42) {

  ## This function conducts log spline maximum
  ## likelihood. Multistarting is supported as is breaking out to
  ## potentially avoid wasted computation (be careful when using this,
  ## however, as it is prone to stopping early).

  ## Save seed prior to setting
  
  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }
  
  set.seed(random.seed)

  penalty <- match.arg(penalty)

  if(missing(x)) stop(" You must provide data")

  ## We set some initial parameters that are placeholders to get
  ## things rolling.

  d.opt <- 0
  s.opt <- 0
  d <- 1
  par.opt <- Inf
  value.opt <- -Inf

  ## Loop through all degrees for every segment starting at
  ## segments.min.

  for(s in segments.min:segments.max) {

    if(do.break & s > 1 & d.opt < d & s.opt < s) break

    ## For smooth densities one can simply restrict degree to at least
    ## 2 (or 3 to be consistent with cubic splines)

    for(d in degree.min:degree.max) {
      if(verbose) cat("\n")
      cat("\r                                                                                                  ")
      cat("\rOptimizing, degree = ",d,", segments = ",s,", degree.opt = ",d.opt, ", segments.opt = ",s.opt," ",sep="")
      dim.p <- dim.bs(basis=basis,degree=d,segments=s)

      ## Multistart if desired.

      for(n in 1:nmulti) {

        ## Can restart to see if we can improve on min... note initial
        ## values totally ad-hoc...

        par.init.out <- par.init(d,s)
        par.init <- par.init.out$par.init
        par.upper <- par.init.out$par.upper

        ## Trap non-convergence, restart from different initial
        ## points, display message if needed (trace>0 up to 6 provides
        ## ever more detailed information for L-BFGS-B)

        optim.out <- list()
        optim.out[[4]] <- 9999
        optim.out$value <- -Inf

        m.attempts <- 0

        while(tryCatch(suppressWarnings(optim.out <- optim(par=par.init,
                                                           fn=sum.log.density,
                                                           gr=if(do.gradient){sum.log.density.gradient}else{NULL},
                                                           upper=par.upper,
                                                           method=method,
                                                           x=x,
                                                           degree=d,
                                                           segments=s,
                                                           er=er,
                                                           penalty=penalty,
                                                           basis=basis,
                                                           knots=knots,
                                                           n.integrate=n.integrate,
                                                           monotone=monotone,
                                                           lbound=lbound,
                                                           ubound=ubound,
                                                           control=list(fnscale=-1,maxit=maxit,if(verbose){trace=1}else{trace=0}))),
                       error = function(e){return(optim.out)})[[4]]!=0 && m.attempts < max.attempts){

          ## If optim fails to converge, reset initial parameters and
          ## try again.

          if(verbose && optim.out[[4]]!=0) {
            if(!is.null(optim.out$message)) cat("\n optim message = ",optim.out$message,sep="")
            cat("\n optim failed (degree = ",d,", segments = ",s,", convergence = ", optim.out[[4]],") re-running with new initial values",sep="")
          }

          par.init.out <- par.init(d,s)
          par.init <- par.init.out$par.init
          par.upper <- par.init.out$par.upper

          m.attempts <- m.attempts+1

        }

        ## Check for a new optimum, overwrite existing values with
        ## new values.

        if(optim.out$value > value.opt) {
          if(verbose && n==1) cat("\n optim improved: d = ",d,", s = ",s,", old = ",formatC(value.opt,format="g",digits=6),", new = ",formatC(optim.out$value,format="g",digits=6),", diff = ",formatC(optim.out$value-value.opt,format="g",digits=6),sep="")

          if(verbose && n>1) cat("\n optim improved (ms ",n,"/",nmulti,"): d = ",d,", s = ",s,", old = ",formatC(value.opt,format="g",digits=6),", new = ",formatC(optim.out$value,format="g",digits=6),", diff = ",formatC(optim.out$value-value.opt,format="g",digits=6),sep="")

          par.opt <- optim.out$par
          d.opt <- d
          s.opt <- s
          value.opt <- optim.out$value
        }

      }

      if(do.break & d.opt < d & s.opt < s) break

    }

  }

  cat("\r                                                                            ")
  if(d.opt==degree.max) warning(paste(" optimal degree equals search maximum (", degree.max,"): rerun with larger degree.max",sep=""))
  if(s.opt==segments.max) warning(paste(" optimal segment equals search maximum (", segments.max,"): rerun with larger segments.max",sep=""))
  if(par.opt[1]>0|par.opt[length(par.opt)]>0) warning(" optim() delivered a positive weight for linear segment (supposed to be negative)")

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
  return(list(degree=d.opt,segments=s.opt,beta=par.opt,fv=value.opt))

}

summary.clsd <- function(object,
                         ...) {

  cat("\nCategorical Logspline Density\n",sep="")
  cat(paste("\nModel penalty: ", format(object$penalty), sep=""))
  cat(paste("\nModel degree/segments: ", format(object$degree),"/",format(object$segments), sep=""))
  cat(paste("\nKnot type: ", format(object$knots), sep=""))
  cat(paste("\nBasis type: ",format(object$basis),sep=""))
  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  cat(paste("\nLog-likelihood: ", format(object$logl), sep=""))
  cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))

  cat("\n\n")

}

print.clsd <- function(x,...)
{
        summary.clsd(x)
}

plot.clsd <- function(object,
                      er=FALSE,
                      distribution=FALSE,
                      ylim,
                      ...) {

  if(!er) {
    order.x <- order(object$x)
    if(distribution){y <- object$distribution[order.x]}else{y <- object$density[order.x]}
    x <- plot(object$x[order.x],
              y,
              ylim=if(missing(ylim)){c(0,max(y))}else{ylim},
              ylab=if(distribution){"Distribution"}else{"Density"},
              xlab="Data",
              type="l",
              ...)
  } else {
    order.xer <- order(object$xer)
    if(distribution){y <- object$distribution.er[order.xer]}else{y <- object$density.er[order.xer]}
    xer <- plot(object$xer[order.xer],
                  y,
                  ylim=if(missing(ylim)){c(0,max(y))}else{ylim},
                  ylab=if(distribution){"Distribution"}else{"Density"},
                  xlab="Data",
                  type="l",
                  ...)
  }
  

}
