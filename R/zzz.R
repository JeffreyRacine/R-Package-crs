.onAttach <- function (lib, pkg) {
  cat("Categorical Regression Splines (version 0.13-7)\n");
  if(is.null(options('crs.messages')$crs.messages))
    options(crs.messages = TRUE)

}
.Last.lib <- function (lpath){
  library.dynam.unload("crs", libpath=lpath) 
  # cat("np unloaded\n")
}
