## This illustration considers the `radial function' and plots the
## results by constructing a 3D real-time rendering plot using OpenGL.

require(crs)
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = TRUE)
require(rgl)

set.seed(42)

.crs_demo_numeric <- function(prompt, default, env = NULL, choices = NULL) {
  value <- if(!is.null(env)) Sys.getenv(env, unset = "") else ""
  if(!nzchar(value) && interactive()) value <- readline(prompt = prompt)
  if(!nzchar(value)) value <- as.character(default)
  value <- suppressWarnings(as.numeric(value))
  if(!is.finite(value) || (!is.null(choices) && !(value %in% choices))) {
    value <- default
  }
  value
}

## Interactively request number of observations, whether to do NOMAD
## or exhaustive search, and if NOMAD the number of multistarts

n <- .crs_demo_numeric("Input the number of observations desired: ",
                       100,
                       "CRS_DEMO_N")
cv <- .crs_demo_numeric("Input the cv method (0=nomad, 1=exhaustive): ",
                        0,
                        "CRS_DEMO_CV",
                        choices = c(0, 1))
cv <- ifelse(cv==0,"nomad","exhaustive")
nmulti <- 1
if(cv=="nomad") {
  nmulti <- .crs_demo_numeric("Input the number of multistarts desired (e.g. 10): ",
                              1,
                              "CRS_DEMO_NMULTI")
}

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

model <- crs(y~x1+x2,
             cv=cv,
             complexity="degree-knots",
             knots="uniform",
             cv.func="cv.aic",
             nmulti=nmulti)

summary(model)

## Create an interactive rgl perspective plot.

plot(model,perspective=TRUE,renderer="rgl",data_rug=TRUE)

## You could animate the results for 15 seconds using the line
## play3d(spin3d(axis=c(0,0,1), rpm=5), duration=15)
## By default you can manually rotate the figure by dragging the plot
## via your mouse/keypad

## Note - to save an rgl figure first get it oriented how you want
## (i.e. resize, rotate etc.) and then call rgl.postscript("foo.pdf","pdf")
## or rgl.snapshot("foo.png").

## For Quarto or R Markdown, place the call in a regular R chunk such
## as the following illustration:
## ```{r}
## x <- rnorm(100); y <- rnorm(100); z <- rnorm(100)
## plot3d(x, y, z)
## ```
