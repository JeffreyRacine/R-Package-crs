.onAttach <- function (lib, pkg) {
	packageStartupMessage("Categorical Regression Splines (version 0.15-5)", domain = NULL,  appendLF = TRUE)
  if(is.null(options('crs.messages')$crs.messages))
    options(crs.messages = TRUE)

}
.Last.lib <- function (lpath){
  library.dynam.unload("crs", libpath=lpath) 
  # cat("np unloaded\n")
}
