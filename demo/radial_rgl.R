require(crs)
require(rgl)

set.seed(42)

## Interactively request number of observations, whether to do NOMAD
## or exhaustive search, and if NOMAD the number of multistarts

n <- as.numeric(readline(prompt="Input the number of observations desired: "))
cv <- as.numeric(readline(prompt="Input the cv method (0=nomad, 1=exhaustive): "))
cv <- ifelse(cv==0,"nomad","exhaustive")
if(cv=="nomad") nmulti <- as.numeric(readline(prompt="Input the number of multistarts desired (e.g. 10): "))
num.eval <- as.numeric(readline(prompt="Input the number of evaluation observations desired (e.g. 50): "))
relief <- as.numeric(readline(prompt="Input the relief exaggeration (e.g. 5): "))

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

model <- crs(y~x1+x2,
             basis="auto",
             cv=cv,
             complexity="degree-knots",
             knots="uniform",
             deriv=1,
             cv.func="cv.aic",
             nmulti=nmulti)

summary(model)

# Perspective plot via rgl (need to also assign colors)
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])
z <- matrix(predict(model,newdata=newdata),num.eval,num.eval)
## For rgl you may want to amplify/reduce the relief (z-magnification)
## and choose a color palette
z <- relief*z
zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- topo.colors(zlen) 
col <- colorlut[ z-zlim[1]+1 ]
## Open an rgl 3d window and use `persp3d', a high-level function for
## 3D surfaces
open3d()
persp3d(x=x1.seq,y=x2.seq,z=z,
        xlab="x1",ylab="x2",zlab="y",
        ticktype="detailed",      
        border="red",
        color=col,
        back="lines",
        main="Conditional Mean")
## Animate the results spinning for 15 seconds... you can manually
## rotate the figure by dragging the plot via your mouse/keypad
## play3d(spin3d(axis=c(0,0,1), rpm=5), duration=15)
