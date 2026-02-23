#
# File:   get.option.types.R
# Author: Zhenghua Nie
# Date:   Mon 16 May 2011
#
# We use ipoptr developed by Jelmer Ypma as the prototype of this package.
# Some code is copied and edited from ipoptr.
# Please reference the license of ipoptr.
#
# This function converts a list with nomad options into
# three sub-lists, where the options are sorted into
# the different value types (integer, numeric, string).
#
# Input: list of nomad options and their values
# Output: list containing three sub-lists by type with nomad options and their values
#
# Copyright (C) 2011 Zhenghua Nie. All Rights Reserved.
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

get.option.types <- function(opts) {

  # initialize list with options sorted by type
  converted.opts <- list( "integer"=list(), "string"=list(), "numeric"=list() )

  # check if we have at least 1 element in the list, otherwise the
  # loop runs from 1 to down 0 and we get errors
  if ( length( opts ) > 0 ) {

    # loop over all options and give them the correct type
    for ( i in 1:length( opts ) ) {
      # Keep all options as strings and let NOMAD parse values.
      # This supports the full NOMAD4 option set without maintaining a hard-coded list.
      tmp.type <- "string"

      if ( tmp.type=="string" ) {
        converted.opts$string[[ names(opts)[i] ]] <- as.character(opts[[i]])
      } else if ( tmp.type=="integer" ) {
        converted.opts$integer[[ names(opts)[i] ]] <- as.integer(opts[[i]])
      } else if ( tmp.type=="numeric" ) {
        converted.opts$numeric[[ names(opts)[i] ]] <- as.numeric(opts[[i]])
      } else {
        stop(paste("Type of option ", names(opts)[i], " not recognized"))
      }
    }
  }

  return ( converted.opts )
}
