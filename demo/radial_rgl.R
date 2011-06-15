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
relief <- as.numeric(readline(prompt="Input the relief exxageration (e.g. 5): "))

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

# Perspective plot
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])
z <- matrix(predict(model,newdata=newdata),num.eval,num.eval)

z <- relief*z ##  # Exaggerate the relief
zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1
#colorlut <- terrain.colors(zlen) # height color lookup table
colorlut <- heat.colors(zlen) # height color lookup table
col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
open3d()
surface3d(x1.seq, x2.seq, z, color=col, back="lines")
