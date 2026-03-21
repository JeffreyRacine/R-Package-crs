.onAttach <- function (lib, pkg) {
  packageStartupMessage(
    sprintf(
      'crs %s: vignette("crs_getting_started", package = "crs")',
      utils::packageDescription(pkg, fields = "Version")
    ),
    domain = NULL,
    appendLF = TRUE
  )
}

.onLoad <- function (lib, pkg) {
  if(is.null(getOption("crs.messages")))
    options(crs.messages = TRUE)
  ## The package np is declared in NAMESPACE and Imports
  ## (DESCRIPTION) to support npglpreg() which eventually will
  ## migrate to the np package, at which time the following are no
  ## longer required.
  if(is.null(getOption("np.messages")))
    options(np.messages = TRUE)
  if(is.null(getOption("np.tree")))
    options(np.tree = FALSE)
}

.onUnload <- function (lpath){
  library.dynam.unload("crs", libpath=lpath)
}
