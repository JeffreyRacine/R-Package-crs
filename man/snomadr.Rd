\name{snomadr}
\alias{snomadr}
\title{
R interface to NOMAD
}
\description{
\code{snomadr} is an R interface to NOMAD (Nonsmooth Optimization by
Mesh Adaptive Direct Search, Abramson, Audet, Couture and Le Digabel
(2011)), an open source software C++ implementation of the Mesh Adaptive
Direct Search (MADS, Le Digabel (2011)) algorithm designed for
constrained optimization of blackbox functions.

NOMAD is designed to find (local) solutions of mathematical optimization
problems of the form

  \preformatted{
   min     f(x)
x in R^n

s.t.       g(x) <= 0 
           x_L <=  x   <= x_U
}           
           
where \eqn{f(x)\colon R^n \to R^k}{f(x): R^n --> R^k} is the objective
function, and \eqn{g(x)\colon R^n \to R^m}{g(x): R^n --> R^m} are
the constraint functions. The vectors \eqn{x_L}{x_L} and \eqn{x_U}{x_U}
are the bounds on the variables \eqn{x}. The functions \eqn{f(x)} and
\eqn{g(x)} can be nonlinear and nonconvex. The variables can be integer,
continuous real number, binary, and categorical.

Kindly see \url{http://www.gerad.ca/nomad/Project/Home.html} and the
references below for details.

}

\usage{
snomadr(eval.f, 
        n, 
        bbin = NULL,  
        bbout = NULL, 
        x0 = NULL, 
        lb = NULL, 
        ub = NULL, 
        nmulti = 0, 
        opts = list(),
        print.output = TRUE, 
        information = list(), 
        snomadr.environment = new.env(), 
        ... ) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{eval.f}{
        function that returns the value of the objective function.
    }
    \item{n}{ 
        the number of variables.
    }

    \item{bbin}{
        types of variables. 0 (CONTINUOUS),  1 (INTEGER),  2 (CATEGORICAL),
        3 (BINARY).
    }
    
    \item{bbout}{
      types of output of \code{eval.f}.  See the NOMAD User Guide
      \url{http://www.gerad.ca/NOMAD/Downloads/user_guide.pdf}.
    }
    
    \item{x0}{
      vector with starting values for the optimization. If it is
        provided and nmulti is bigger than 1, \code{x0} will be the first
        initial point for multiple initial points.
      }
    \item{lb}{
        vector with lower bounds of the controls (use -1.0e19 for controls without lower bound).
    }
    \item{ub}{
        vector with upper bounds of the controls (use  1.0e19 for controls without upper bound).
    }
    \item{nmulti}{
      when it is missing, or it is equal to 0 and \code{x0} is provided,
        \code{snomadRSolve} will be called to solve the problem. Otherwise,
        \code{smultinomadRSolve} will be called.
      }
    \item{opts}{
      
      list of options for NOMAD, see the NOMAD user guide
      \url{http://www.gerad.ca/NOMAD/Downloads/user_guide.pdf}.  Options
      can also be set by nomad.opt which should be in the folder where R
      (\code{snomadr}) is executed. Options that affect the solution and
      their defaults and some potential values are
      
      \code{"MAX_BB_EVAL"=10000}

      \code{"INITIAL_MESH_SIZE"=1}

      \code{"MIN_MESH_SIZE"="r1.0e-10"}
      
      \code{"MIN_POLL_SIZE"="r1.0e-10"}

      Note that the \code{"r..."} denotes relative measurement (relative
      to \code{degree.min} and \code{degree.max} and \code{segments.min}
      and \code{segments.max}).

      Note that decreasing the maximum number of black box evaluations
      (\code{"MAX_BB_EVAL"}) will terminate search sooner and may result
      in a less accurate solution. For complicated problems you may want
      to increase this value. When experimenting it is desirable to set
      \code{"DISPLAY_DEGREE"=1} (or a larger integer) to get some sense
      for how the algorithm is progressing.

      }

    \item{print.output}{
      when FALSE, no output from \code{snomadr} is displayed on
      the screen. If \code{"DISPLAY_DEGREE"=0,} the option of NOMAD, there will
      be no output from NOMAD. Higher integer values provide
      successively more detail.
    }

		\item{information}{
			is a list. \code{snomadr} will call \code{snomadRInfo} to return
			the information about NOMAD according to the values of \code{"info"},
			\code{"version"} and \code{"help"}.

			\code{"info"="-i"}: display the usage and copyright of NOMAD. 

			\code{"version"="-v"}: display the version of NOMAD you are using.

			\code{"help"="-h"}: display all options. You also can display a specific option, 
			for example, \code{"help"="-h x0"},  this will tell you how to set \code{x0}.
		}

    \item{snomadr.environment}{
        environment that is used to evaluate the functions. Use this to pass 
        additional data or parameters to a function.
    }
    \item{...}{
      arguments that will be passed to the user-defined objective and
      constraints functions. See the examples below.
    }
  }
  \details{
    \code{snomadr} is used in the \pkg{crs} package to numerically
    minimize an objective function with respect to the spline degree,
    number of knots, and optionally the kernel bandwidths when using
    \code{\link{crs}} with the option \code{cv="nomad"} (default). See
    \code{\link{frscvNOMAD}} and \code{\link{krscvNOMAD}} for the
    functions called when \code{cv="nomad"} while using
    \code{\link{crs}}.

    However, note that for simple problems involving one predictor
    exhaustive search may be faster and more accurate, so bear in mind
    that \code{cv="exhaustive"} can be useful when using
    \code{\link{crs}}. Naturally, exhaustive search is also useful for
    verifying solutions returned by \code{snomadr}. See
    \code{\link{frscv}} and \code{\link{krscv}} for the functions called
    when \code{cv="exhaustive"} while using \code{\link{crs}}.
  }
\value{
    The return value contains a list with the inputs, and additional elements
    \item{call}{the call that was made to solve}
    \item{status}{integer value with the status of the optimization }
    \item{message}{more informative message with the status of the optimization}
    \item{iterations}{number of iterations that were executed, if multiple initial points are set, 
        this number will be the sum for each initial point.}
    \item{objective}{value if the objective function in the solution}
    \item{solution}{optimal value of the controls}
}
\references{
  M.A. Abramson, C. Audet, G. Couture, J.E. Dennis, Jr., and S. Le
  Digabel (2011), \dQuote{The NOMAD project}. Software available at
  http://www.gerad.ca/nomad.

  S. Le Digabel (2011), \dQuote{Algorithm 909: NOMAD: Nonlinear
  optimization with the MADS algorithm}. ACM Transactions on
  Mathematical Software, 37(4):44:1-44:15.
  
}
\author{
  Zhenghua Nie <niez@mcmaster.ca>
}
\seealso{
    \code{\link{optim}}
    \code{\link{nlm}}
    \code{\link{nlminb}}
}
\examples{
# list all options
snomadr(information=list("help"="-h"))

# print given option,  for example,  MESH_SIZE
snomadr(information=list("help"="-h MESH_SIZE"))

#print the version of NOMAD
snomadr(information=list("version"="-v"))

#print usage and copyright
snomadr(information=list("info"="-i"))

# This is the example found in
#  NOMAD/examples/basic/library/single_obj/basic_lib.cpp

eval.f <- function ( x ) {

    f <- c(Inf, Inf, Inf);
    n <- length (x);

    if ( n == 5 && ( is.double(x) || is.integer(x) ) ) {
        f[1] <- x[5];
        f[2] <- sum ( (x-1)^2 ) - 25;
        f[3] <- 25 - sum ( (x+1)^2 );
    }  

    return ( as.double(f) );
}

# initial value
x0 <- rep(0.0, 5 )

bbin <-c(1, 1, 1, 1, 1)
#bounds
lb <- rep(-6.0,5 )
ub <- c(5.0, 6.0, 7.0, 1000000, 100000)

bbout <-c(0, 2, 1)
#options
opts <-list("MAX_BB_EVAL"=500,
            "MIN_MESH_SIZE"=0.001,
            "INITIAL_MESH_SIZE"=0.1,
            "MIN_POLL_SIZE"=0.0001)

snomadr(eval.f=eval.f,n=5,  x0=x0, bbin=bbin, bbout=bbout, lb=lb, ub=ub, opts=opts)


# how to transfer other parameters into eval.f
#
# First example: supply additional arguments in user-defined functions
#

# objective function and gradient in terms of parameters
eval.f.ex1 <- function(x, params) { 
    return( params[1]*x^2 + params[2]*x + params[3] ) 
}

# define parameters that we want to use
params <- c(1,2,3)

# define initial value of the optimization problem
x0 <- 0

# solve using snomadr 
snomadr( n          =1, 
        x0          = x0, 
        eval.f      = eval.f.ex1, 
        params      = params )


#
# Second example: define an environment that contains extra parameters
#

# objective function and gradient in terms of parameters
# without supplying params as an argument
eval.f.ex2 <- function(x) { 
    return( params[1]*x^2 + params[2]*x + params[3] ) 
}

# define initial value of the optimization problem
x0 <- 0

# define a new environment that contains params
auxdata        <- new.env()
auxdata$params <- c(1,2,3)

# pass the environment that should be used to evaluate functions to snomadr 
snomadr( n                 =1, 
        x0                 = x0, 
        eval.f             = eval.f.ex2, 
        snomadr.environment = auxdata )

# solve using algebra
cat( paste( "Minimizing f(x) = ax^2 + bx + c\n" ) )
cat( paste( "Optimal value of control is -b/(2a) = ", -params[2]/(2*params[1]), "\n" ) )
cat( paste( "With value of the objective function f(-b/(2a)) = ", eval.f.ex1( -params[2]/(2*params[1]), params ), "\n" ) )


# The following example is NOMAD/examples/advanced/multi_start/multi.cpp
# This will call smultinomadRSolve to resolve the problem.  
eval.f.ex1 <- function(x, params) { 
    M<-as.numeric(params$M)
    NBC<-as.numeric(params$NBC)

    f<-rep(0, M+1)
    x<-as.numeric(x)

    x1 <- rep(0.0, NBC)
    y1 <- rep(0.0, NBC)

    x1[1]<-x[1]
    x1[2]<-x[2]
    y1[3]<-x[3]
    x1[4]<-x[4]
    y1[4]<-x[5]

    epi <- 6
    for(i in 5:NBC){
        x1[i]<-x[epi]
        epi <- epi+1
        y1[i]<-x[epi]
        epi<-epi+1
    }
    constraint <- 0.0
    ic <- 1
    f[ic]<-constraint
    ic <- ic+1

    constraint <- as.numeric(1.0)
    distmax <- as.numeric(0.0)
    avg_dist <- as.numeric(0.0)
    dist1<-as.numeric(0.0)

    for(i in 1:(NBC-1)){
        for (j in (i+1):NBC){
            dist1 <- as.numeric((x1[i]-x1[j])*(x1[i]-x1[j])+(y1[i]-y1[j])*(y1[i]-y1[j]))
            
            if((dist1 > distmax)) {distmax <- as.numeric(dist1)}
            if((dist1[1]) < 1) {constraint <- constraint*sqrt(dist1)}
            else if((dist1) > 14) {avg_dist <- avg_dist+sqrt(dist1)}
        }
    }

    if(constraint < 0.9999) constraint <- 1001.0-constraint
    else constraint = sqrt(distmax)+avg_dist/(10.0*NBC)

    f[2] <- 0.0
    f[M+1] <- constraint 


    return(as.numeric(f) ) 
}

# define parameters that we want to use
params<-list()
NBC <- 5
M <- 2
n<-2*NBC-3

params$NBC<-NBC
params$M<-M
x0<-rep(0.1, n)
lb<-rep(0, n)
ub<-rep(4.5, n)

eval.f.ex1(x0, params)

bbout<-c(2, 2, 0)
nmulti=5
bbin<-rep(0, n)
# define initial value of the optimization problem

# solve using snomadRSolve
snomadr(n            = as.integer(n), 
        x0           = x0, 
        eval.f       = eval.f.ex1, 
        bbin         = bbin, 
        bbout        = bbout, 
        lb           = lb, 
        ub           = ub, 
        params       = params )

# solve using smultinomadRSolve, if x0 is provided,  x0 will
# be the first initial point,  otherwise,  the program will
# check best_x.txt,  if it exists,  it will be read in as
# the first initial point. Other initial points will be
# generated by uniform distribution.
# nmulti represents the number of mads will run.
#
snomadr(n            = as.integer(n), 
        eval.f       = eval.f.ex1, 
        bbin         = bbin, 
        bbout        = bbout, 
        lb           = lb, 
        ub           = ub, 
        nmulti = as.integer(nmulti), 
        print.output = TRUE, 
        params       = params )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
\keyword{ interface }

