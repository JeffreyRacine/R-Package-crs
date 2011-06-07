# Copyright (C) 2011 Zhenghua  Nie. All Rights Reserved.
# This code is published under the Common Public License.
#
# File:   parameters.R
# Author: Zhenghua Nie
# Date:   Mon 16 May 2011
#
# We use ipoptr developed by Jelmer Ypma as the prototype of this package.
# Some code is copied and edited from ipoptr. 
# Please reference the license of ipoptr.
#
# Example showing two ways how we can have an objective 
# function depend on parameters or data. The objective
# function is a simple polynomial.

require('crs')


#
# First example: supply additional arguments in user-defined functions
#

# objective function and gradient in terms of parameters
eval_f_ex1 <- function(x, params) { 
    return( params[1]*x^2 + params[2]*x + params[3] ) 
}

# define parameters that we want to use
params <- c(1,2,3)

# define initial value of the optimzation problem
x0 <- 0

# solve using ipoptr
snomadr( n          =1, 
				x0          = x0, 
        eval_f      = eval_f_ex1, 
        params      = params )


#
# Second example: define an environment that contains extra parameters
#

# objective function and gradient in terms of parameters
# without supplying params as an argument
eval_f_ex2 <- function(x) { 
    return( params[1]*x^2 + params[2]*x + params[3] ) 
}

# define initial value of the optimzation problem
x0 <- 0

# define a new environment that contains params
auxdata        <- new.env()
auxdata$params <- c(1,2,3)

# pass the environment that should be used to evaluate functions to ipoptr
snomadr( n                 =1, 
				x0                 = x0, 
        eval_f             = eval_f_ex2, 
        snomadr_environment = auxdata )


# solve using algebra
cat( paste( "Minimizing f(x) = ax^2 + bx + c\n" ) )
cat( paste( "Optimal value of control is -b/(2a) = ", -params[2]/(2*params[1]), "\n" ) )
cat( paste( "With value of the objective function f(-b/(2a)) = ", eval_f_ex1( -params[2]/(2*params[1]), params ), "\n" ) )
