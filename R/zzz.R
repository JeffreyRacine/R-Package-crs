.onAttach <- function (lib, pkg) {
	packageStartupMessage("Categorical Regression Splines (version 0.15-20)\n[vignette(\"crs_faq\") provides answers to frequently asked questions]", domain = NULL,  appendLF = TRUE)
  if(is.null(options('crs.messages')$crs.messages))
    options(crs.messages = TRUE)

}
.onUnload <- function (lpath){
  library.dynam.unload("crs", libpath=lpath)
}
