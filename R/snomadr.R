# File:   snomadr.R
# Author: Zhenghua Nie
# Date:   Mon 16 May 2011
#
# We use ipoptr developed by Jelmer Ypma as the prototype of this package.
# Some code is copied and edited from ipoptr.
# Please reference the license of ipoptr.
#
# Input:
#   n  : the number of variables
#    x0 : vector with initial values
#    eval.f : function to evaluate objective function
#    lb : lower bounds of the control
#    ub : upper bounds of the control
#    opts : list with options that are passed to NOMAD
#       ... : arguments that will be passed to user-defined functions
#
# Output: structure with inputs and
#    call : the call that was made to solve
#    status : integer value with the status of the optimization (0 is success)
#    message : more informative message with the status of the optimization
#    bbe : number of the objective function that were executed
#    iterations : number of iterations that were executed
#    objective : value if the objective function in the solution
#    solution : optimal value of the controls
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

nomad4.mads.defaults <- function() {
  # Default MADS profile using NOMAD4 option names only.
  # User-supplied opts always take precedence.
  list(
    "QUAD_MODEL_SEARCH" = "yes",
    "SGTELIB_MODEL_SEARCH" = "no",
    "NM_SEARCH" = "yes",
    "SPECULATIVE_SEARCH" = "yes",
    "EVAL_OPPORTUNISTIC" = "yes",
    "EVAL_QUEUE_SORT" = "QUADRATIC_MODEL",
    "DIRECTION_TYPE" = "ORTHO N+1 QUAD",
    # Tuned default for quad-model search box factor.
    "QUAD_MODEL_SEARCH_BOX_FACTOR" = "2.0",
    "QUAD_MODEL_BOX_FACTOR" = "2.0"
  )
}

merge.nomad4.mads.defaults <- function(opts) {
  defaults <- nomad4.mads.defaults()
  if (length(opts) == 0) {
    return(defaults)
  }
  for (k in names(defaults)) {
    if (is.null(opts[[k]])) {
      opts[[k]] <- defaults[[k]]
    }
  }
  opts
}

nomad4.fr.defaults <- function() {
  # Path-specific defaults for frscvNOMAD.
  # Empirically reduces runtime while preserving objective parity.
  list(
    "QUAD_MODEL_SEARCH" = "no",
    "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS",
    "SIMPLE_LINE_SEARCH" = "yes",
    "SPECULATIVE_SEARCH" = "no",
    "DIRECTION_TYPE" = "ORTHO N+1 NEG"
  )
}

merge.nomad4.fr.defaults <- function(opts) {
  defaults <- nomad4.fr.defaults()
  if (length(opts) == 0) {
    return(defaults)
  }
  for (k in names(defaults)) {
    if (is.null(opts[[k]])) {
      opts[[k]] <- defaults[[k]]
    }
  }
  opts
}

nomad4.kr.defaults <- function() {
  # Path-specific defaults for krscvNOMAD.
  # Aggressive speed-oriented defaults.
  list(
    "QUAD_MODEL_SEARCH" = "no",
    "EVAL_QUEUE_SORT" = "DIR_LAST_SUCCESS",
    "DIRECTION_TYPE" = "ORTHO 2N"
  )
}

merge.nomad4.kr.defaults <- function(opts) {
  defaults <- nomad4.kr.defaults()
  if (length(opts) == 0) {
    return(defaults)
  }
  for (k in names(defaults)) {
    if (is.null(opts[[k]])) {
      opts[[k]] <- defaults[[k]]
    }
  }
  opts
}

.crs_nomad_option_value <- function(opts, key) {
  opt.names <- names(opts)
  if (is.null(opt.names))
    return(NULL)
  key <- toupper(key)
  opt.index <- which(toupper(opt.names) == key)
  if (!length(opt.index))
    return(NULL)
  opts[[opt.index[1L]]]
}

.crs_nomad_has_option <- function(opts, key) {
  !is.null(.crs_nomad_option_value(opts, key))
}

.crs_nomad_option_text <- function(x) {
  if (is.null(x)) return(NULL)
  paste(as.character(x), collapse = " ")
}

.crs_nomad_apply_eval_budget <- function(opts,
                                         max.bb.eval = NULL,
                                         max.eval = NULL,
                                         default.max.eval = NULL,
                                         context = "NOMAD") {
  if (!.crs_nomad_has_option(opts, "MAX_BB_EVAL") && !is.null(max.bb.eval))
    opts$"MAX_BB_EVAL" <- max.bb.eval

  if (!is.null(max.eval)) {
    old <- .crs_nomad_option_value(opts, "MAX_EVAL")
    if (!is.null(old) &&
        !identical(.crs_nomad_option_text(old),
                   .crs_nomad_option_text(max.eval))) {
      stop(paste0(
        context,
        " received conflicting MAX_EVAL controls; use either max.eval or ",
        "opts$MAX_EVAL, not both with different values"
      ), call. = FALSE)
    }
    opts$"MAX_EVAL" <- max.eval
  } else if (!.crs_nomad_has_option(opts, "MAX_EVAL") &&
             !is.null(default.max.eval)) {
    opts$"MAX_EVAL" <- default.max.eval
  } else if (!.crs_nomad_has_option(opts, "MAX_EVAL")) {
    max.bb <- .crs_nomad_option_value(opts, "MAX_BB_EVAL")
    max.bb.num <- suppressWarnings(as.numeric(max.bb))
    if (length(max.bb.num) == 1L && is.finite(max.bb.num) && max.bb.num > 0)
      opts$"MAX_EVAL" <- max.bb
  }

  opts
}

.crs_nomad_effective_options <- function(opts) {
  keys <- c("MAX_BB_EVAL",
            "MAX_EVAL",
            "DIRECTION_TYPE",
            "EVAL_QUEUE_SORT",
            "EVAL_OPPORTUNISTIC",
            "SPECULATIVE_SEARCH",
            "SIMPLE_LINE_SEARCH",
            "QUAD_MODEL_SEARCH",
            "SGTELIB_MODEL_SEARCH",
            "NM_SEARCH",
            "INITIAL_MESH_SIZE",
            "MIN_MESH_SIZE",
            "MIN_FRAME_SIZE",
            "INITIAL_FRAME_SIZE")

  out <- list()
  for (key in keys) {
    value <- .crs_nomad_option_value(opts, key)
    if (!is.null(value))
      out[[key]] <- value
  }
  out
}

.crs_nomad_attach_effective_options <- function(summary, opts) {
  if (is.null(summary)) return(NULL)
  summary$effective.options <- .crs_nomad_effective_options(opts)
  summary
}

.crs_nomad_display_degree <- function(opts) {
  value <- .crs_nomad_option_value(opts, "DISPLAY_DEGREE")
  if (is.null(value) || length(value) != 1L)
    return(NULL)
  value <- suppressWarnings(as.integer(value))
  if (is.na(value))
    return(NULL)
  value
}

.crs_nomad_stop_unsafe_display_degree <- function(display.degree) {
  message <- paste(
    "NOMAD option DISPLAY_DEGREE >= 3 is currently unsafe in the embedded",
    "snomadr() route because NOMAD4's native output queue can abort the R",
    "session after repeated solves. Use DISPLAY_DEGREE 0, 1, or 2; keep the",
    "default crs single-line progress display for routine fits.",
    sep = " "
  )
  condition <- structure(
    list(message = message,
         call = NULL,
         display.degree = display.degree),
    class = c("crs_error_nomad_display_degree", "error", "condition")
  )
  stop(condition)
}

.crs_nomad_check_display_degree <- function(opts) {
  display.degree <- .crs_nomad_display_degree(opts)
  if (!is.null(display.degree) && display.degree >= 3L)
    .crs_nomad_stop_unsafe_display_degree(display.degree)
  invisible(opts)
}

snomadr <-
  function( eval.f,
            n,
            bbin = NULL,
            bbout = NULL,
            x0 = NULL,
            lb = NULL,
            ub = NULL,
            nmulti = 0,   #0: call single nomad,
            random.seed = 0,  # seed will be used for generating multiple initial points.
            opts = list(),
            display.nomad.progress = TRUE,  #0: if it is FALSE,  there will be no output in snomadr, if DISPLAY_DEGREE=0 and print_output is true, there will be no any output.
            information = list(),
            snomadr.environment = new.env(),
            ... ) {

    ## Save seed prior to setting

    seed.state <- .crs_capture_seed()

    set.seed(random.seed)

    if(length(information) > 0) {
      sinfo <- NULL
      sversion <- NULL
      shelp <- NULL
      if(!is.null(information$info)) sinfo <-information$info
      if(!is.null(information$version)) sversion <- information$version
      if(!is.null(information$help)) shelp <- information$help

      ret <- list( info=sinfo, version=sversion, help=shelp, snomadr.environment)

      attr(ret, "class") <- "snomadr"

      ## Add the current call to the list
      ret$call <- match.call()

      ## Pass snomadr object to C code
      solution <- .Call( snomadRInfo, ret )
      ## Remove the environment from the return object
      ret$environment <- NULL
      ## We have not implemented the following output from snomadRInfo
      ##ret$Info<-solution$Info
      ##ret$Version<-solution$Version
      ##ret$Help<-solution$Help

      return(ret)

    }

    information <- NULL  ## We will check whether it is NULL in print.snomadr.

    ## The number of variables should not be null.
    if (missing(n ) || missing(eval.f)) stop("Must provide the objective function and the number of variables")
    if(missing(nmulti)||nmulti < 0) nmulti <- 0
    if(missing(display.nomad.progress)) display.nomad.progress <- TRUE
    opts <- merge.nomad4.mads.defaults(opts)
    .crs_nomad_check_display_degree(opts)

    ## Define 'continuous' to types of variables
    if (is.null(bbin) ) { bbin <- rep (0, n)}
    ## Define 'NOMAD::OBJ' to the output type of the function
    if (is.null(bbout)) {bbout <- rep(0, 1)}

    ## Define 'infinite' lower and upper bounds of the control if they haven't been set
    if ( is.null( lb ) ) { lb <- rep( -Inf, n ) }
    if ( is.null( ub ) ) { ub <- rep(  Inf, n ) }

    ## We don't need to generate the initial point for multiple mads runs.
    if(is.null(x0)&&nmulti < 1){
      x0<-rep(0.0, n)
      for(i in seq_len(n)){
        x0[i] <- runif(1, min=lb[i], max=ub[i])
      }
    }

    ## Change the environment of the functions that we're calling the
    ## environment of the eval.f is changed below (if it exists)
    environment( eval.f ) <- snomadr.environment

    ## Internal function to check the arguments of the functions
    checkFunctionArguments <- function( fun, arglist, funname ) {
      if (!is.function(fun)) stop(paste(funname, " must be a function\n", sep = ""))

      fargs <- formals(fun)
      if (length(fargs) <= 1L) {
        if (length(arglist) > 0L) {
          argnames.supplied <- names(arglist)
          if (is.null(argnames.supplied)) {
            argnames.supplied <- rep.int("", length(arglist))
          }
          bad.idx <- 1L
          bad.name <- if (nzchar(argnames.supplied[bad.idx])) {
            argnames.supplied[bad.idx]
          } else {
            paste0("argument #", bad.idx)
          }
          stop(paste("'", bad.name, "' passed to (...) in 'snomadr' but this is not required in the ", funname, " function.\n", sep = ""))
        }
        return(0)
      }

      is_missing_default <- function(x) identical(x, quote(expr = ))

      extra_formals <- fargs[-1L]
      extra_names <- names(extra_formals)
      has_dots <- any(extra_names == "...")
      formal_names <- extra_names[extra_names != "..."]

      argnames.supplied <- names(arglist)
      if (is.null(argnames.supplied)) {
        argnames.supplied <- rep.int("", length(arglist))
      }

      named_idx <- which(nzchar(argnames.supplied))
      unnamed_idx <- which(!nzchar(argnames.supplied))
      named_args <- argnames.supplied[named_idx]

      supplied_by_name <- rep.int(FALSE, length(formal_names))
      names(supplied_by_name) <- formal_names

      if (length(named_args) > 0L) {
        bad_named <- named_args[is.na(match(named_args, formal_names))]
        if (!has_dots && length(bad_named) > 0L) {
          stop(paste("'", bad_named[1L], "' passed to (...) in 'snomadr' but this is not required in the ", funname, " function.\n", sep = ""))
        }
        valid_named <- named_args[!is.na(match(named_args, formal_names))]
        if (length(valid_named) > 0L) {
          supplied_by_name[unique(valid_named)] <- TRUE
        }
      }

      unmatched_formals <- formal_names[!supplied_by_name]
      n_unnamed <- length(unnamed_idx)
      n_positional_match <- min(length(unmatched_formals), n_unnamed)
      supplied_positional <- character(0L)
      if (n_positional_match > 0L) {
        supplied_positional <- unmatched_formals[seq_len(n_positional_match)]
      }

      if (!has_dots && n_unnamed > n_positional_match) {
        bad_pos <- unnamed_idx[n_positional_match + 1L]
        stop(paste("'", paste0("argument #", bad_pos), "' passed to (...) in 'snomadr' but this is not required in the ", funname, " function.\n", sep = ""))
      }

      supplied_all <- rep.int(FALSE, length(formal_names))
      names(supplied_all) <- formal_names
      if (any(supplied_by_name)) {
        supplied_all[names(supplied_by_name)[supplied_by_name]] <- TRUE
      }
      if (length(supplied_positional) > 0L) {
        supplied_all[supplied_positional] <- TRUE
      }

      required_formals <- extra_names[
        extra_names != "..." & vapply(extra_formals[extra_names != "..."], is_missing_default, logical(1L))
      ]
      missing_required <- required_formals[!supplied_all[required_formals]]
      if (length(missing_required) > 0L) {
        stop(paste(funname, " requires argument '", missing_required[1L], "' but this has not been passed to the 'snomadr' function.\n", sep = ""))
      }

      return(0)
    }

    ## Extract list of additional arguments and check user-defined
    ## functions There is an error when building the package crs, so
    ## Zhenghua commented the the following two lines.

    arglist <- list(...)
    checkFunctionArguments( eval.f, arglist, 'eval.f' )

    ## Write wrappers around user-defined functions to pass additional
    ## arguments
    eval.f.wrapper <- function(x){ eval.f(x,...) }

    ## Build snomadr object
    ret <- list("eval.f"=eval.f.wrapper,
                "n"=as.integer(n),
                "bbin"=as.integer(bbin),
                "bbout"=as.integer(bbout),
                "x0"=x0,
                "lower.bounds"=lb,
                "upper.bounds"=ub,
                "nmulti"=as.integer(nmulti),
                "random.seed"=as.integer(random.seed),
                "options"=get.option.types(opts),
                "print.output"=display.nomad.progress,
                "snomadr.environment"=snomadr.environment)

    attr(ret, "class") <- "snomadr"

    ## Add the current call to the list
    ret$call <- match.call()

    ## Pass snomadr object to C code
    if(nmulti == 0){
      solution <- .Call( snomadRSolve, ret )
    } else {
      solution <- .Call( smultinomadRSolve, ret )
    }

    ## Remove the environment from the return object
    ret$environment <- NULL

    ## Add solution variables to object
    ret$status <- solution$status
    ret$message <- solution$message
    ret$bbe <- solution$bbe
    ret$cache.hits <- solution$cache.hits
    ret$cache.size <- solution$cache.size
    ret$callback.evaluations <- solution$callback.evaluations
    ret$total.evaluations <- solution$total.evaluations
    ret$iterations <- solution$iterations
    ret$objective <- solution$objective
    ret$solution <- solution$solution

    ## Restore seed

    .crs_restore_seed(seed.state)

    return( ret )

  }
