.onAttach <- function (lib, pkg) {
  packageStartupMessage(
    sprintf(
      'crs %s: examples and guides at https://jeffreyracine.github.io/gallery/; vignette("crs_getting_started", package = "crs")',
      utils::packageDescription(pkg, fields = "Version")
    ),
    domain = NULL,
    appendLF = TRUE
  )
}

.onLoad <- function (lib, pkg) {
  if(is.null(getOption("crs.messages")))
    options(crs.messages = TRUE)
}

.onUnload <- function (lpath){
  library.dynam.unload("crs", libpath=lpath)
}
