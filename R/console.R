#### single line console management
newLineConsole <- function(parentConsole = NULL){
  newConsole <- list(lineList = list(),
                     numLineElements = 0,
                     tabLen = 8,
                     lineLen = 0,
                     startCol = 1,
                     lineMsg = '',
                     progress = NULL)

  if(!is.null(parentConsole)) {
    newConsole$startCol <- parentConsole$startCol + parentConsole$lineLen
    newConsole$progress <- parentConsole$progress
  }

  if(is.null(newConsole$progress)) {
    newConsole$progress <- .crs_progress_status_begin(surface = "console")
  }
  return(newConsole)
}


charRep <- function(x, times, ...){
  if(times != 0)
    return(rep(x = x, times = times, ...))
  else
    return('')
}

toMsg <- function(msg, console = stop("no console provided")){
  ## attempt to provide minimal support for other data types
  msg <- paste(msg, collapse=' ')

  ## explode the message
  tMsg <- unlist(strsplit(msg, '\\\t'))

  ## add trailing tab if necessary
  if(substr(msg, nchar(msg), nchar(msg)) == '\t')
    tMsg[length(tMsg)+1] = ""

  if(length(tMsg) > 1){
    tmpPos <- console$lineLen + console$startCol

    ## do tabbing
    fills <- character(length(tMsg))
    for(i in seq_along(tMsg)) {
      tmpPos <- tmpPos + nchar(tMsg[i])
      tFill <- paste(charRep(' ', console$tabLen - (tmpPos %% console$tabLen)),
                     collapse = '')
      tmpPos <- tmpPos + nchar(tFill)
      fills[i] <- tFill
    }
    msg <- paste(tMsg, fills, sep = '', collapse = '')
  }
  list(msg = msg,
       len = nchar(msg))
}


printPush <- function(msg, console = stop("no console provided")){
  if(isTRUE(getOption("crs.messages"))){
    console$numLineElements <- console$numLineElements + 1
    console$lineList[console$numLineElements] <- NA
    console$lineList[[console$numLineElements]] <- toMsg(msg = msg, console = console)
    console$lineLen <- console$lineLen + console$lineList[[console$numLineElements]]$len
    console$lineMsg <- paste(console$lineMsg, console$lineList[[console$numLineElements]]$msg, sep = '')
    .crs_progress_status_update(console$progress, console$lineMsg)
  }
  return(console)
}

printPop <- function(console = stop("no console provided")){
  if(console$numLineElements > 0 && isTRUE(getOption("crs.messages"))) {
    console$lineLen <- console$lineLen - console$lineList[[console$numLineElements]]$len
    console$lineMsg <- substr(console$lineMsg, 1, console$lineLen)
    stopifnot(console$lineLen >= 0)
    console$lineList[console$numLineElements] <- NULL
    console$numLineElements <- console$numLineElements - 1
    if(console$numLineElements > 0 && nzchar(console$lineMsg)) {
      .crs_progress_status_update(console$progress, console$lineMsg)
    } else {
      .crs_progress_status_clear(console$progress)
    }
  }
  return(console)
}

printClear <- function(console = stop("no console provided")){
  if(console$numLineElements > 0 && isTRUE(getOption("crs.messages"))){
    .crs_progress_status_clear(console$progress)
    console$lineLen <- 0
    console <- newLineConsole(console)
  }
  return(console)
}
