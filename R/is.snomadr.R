#
# File:   is.snomadr.R
# Author: Zhenghua Nie
# Date:   Mon 16 May 2011
#
# We use ipoptr developed by Jelmer Ypma as the prototype of this package.
# Some code is copied and edited from ipoptr. 
# Please reference the license of ipoptr.
#
# Input: object
# Output: bool telling whether the object is an snomadr or not
#
# Copyright (C) 2011 Zhenghua  Nie. All Rights Reserved.
# This code is published under GNU GENERAL PUBLIC LICENSE.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,  or
# (at your option) any later version.
#      
# This program is distributed WITHOUT ANY WARRANTY. See the
# GNU General Public License for more details.
#           
# If you do not have a copy of the GNU General Public License,  
# write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

is.snomadr <- function(x) {
    
    # Check whether the object exists and is a list
    if( is.null(x) ) { return( FALSE ) }
    if( !is.list(x) ) { return( FALSE ) }

    # Check whether a correct environment was specified
    stopifnot( is.environment(x$snomadr.environment) )
    
    # Check whether the needed functions are supplied
    stopifnot( is.function(x$eval.f) )
    
    # Check whether bounds are defined for all controls
    stopifnot( length( x$upper.bounds ) == length( x$lower.bounds ) )
    
    # Check whether the initial value is within the bounds
		if(length(x$x0) ==x$n){
    stopifnot( all( x$x0 >= x$lower.bounds ) )
    stopifnot( all( x$x0 <= x$upper.bounds ) )
    # Check the length of some return values
    stopifnot( length(x$eval.f( x$x0 ))==length(x$bbout) )
    
    # Check the whether we don't have NA's in initial values
    stopifnot( all(!is.na(x$eval.f( x$x0 ))) )
		}

    return( TRUE )
}
