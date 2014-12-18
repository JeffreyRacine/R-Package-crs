.onAttach <- function (lib, pkg) {
	packageStartupMessage("Categorical Regression Splines (version 0.15-23)\n[vignette(\"crs_faq\") provides answers to frequently asked questions]", domain = NULL,  appendLF = TRUE)
}

.onLoad <- function (lib, pkg) {
    if(is.null(options('crs.messages')$crs.messages))
    options(crs.messages = TRUE)
}

.onUnload <- function (lpath){
  library.dynam.unload("crs", libpath=lpath)
}
