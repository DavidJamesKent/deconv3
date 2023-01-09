library(devtools)

## Install the package. This works on my Linux computer; not sure if it will
## work on Windows or Mac OS
devtools::install_github("https://github.com/DavidJamesKent/deconv3")

## Demo:

## Suppose our unknown density is an awkward mixture of a Gamma and a normal:
fx <- function(x) 0.5*dgamma(x,5,2) + 0.5*dnorm(x,5,2)
rx <- function(n) ifelse(rbinom(n,1,0.5),rgamma(n,5,2),rnorm(n,5,2))
hist(rx(1e5),breaks=200,freq=F,border=F)
curve(fx(x),add=T,col=2,n=1e3)
m <- integrate(function(x) x*fx(x),-Inf,Inf)$value
v <- integrate(function(x) (x - m)^2*fx(x),-Inf,Inf)$value

## Draw a sample of size n = 500 from this distribution and pollute
## with Gaussian errors that give 10% measurement error.
set.seed(1)
n <- 500
X1 <- rx(n)
mev <- v/9 # Measurement error variance
E1 <- rnorm(n,mean=0,sd=sqrt(mev))
Y1 <- X1+E1

## The prep to use the package: we need to specify the characteristic
## function of the measurement error density (called "g" in the
## paper), which can be done with this helper function:
gtwid <- deconv3::getgtwid(g="Gaussian",param=c(s=sqrt(mev)))

## - K is the number of bins in the histogram estimate of the density of Y
## - qdim is the dimension of the spline basis used for representing the
## estimated density
## - If we do not supply a penalty parameter lambda, it uses an iterated
## bootstrap method.
result1 <- deconv3::decon(Y=Y1,gtwid=gtwid,K=100,qdim=50)
## The result is a function stored in `fhat`:
curve(fx(x),-5,15,n=1e3,ylim=c(0,0.3))
curve(result1$fhat(x),add=T,col=2,n=1e3)

## If we plan on repeatedly calling decon() with similar data, we can avoid
## re-computing the miscellaneous bits (inner products, convolved spline basis
## functions, etc.) by passing the `control` object from a previous call:
X2 <- rx(n)
E2 <- rnorm(n,mean=0,sd=sqrt(mev))
Y2 <- X2+E2
result2 <- deconv3::decon(Y=Y2,control=result1$control)
curve(fx(x),-5,15,n=1e3,ylim=c(0,0.3),lwd=2)
curve(result1$fhat(x),add=T,col=2,n=1e3)
curve(result2$fhat(x),add=T,col=4,n=1e3)
legend("topright",legend=c("Unknown Target","Rep 1","Rep 2"),lty=1,lwd=c(2,1,1),col=c(1,2,4))

## We can see what was chosen for the penalty parameters by extracting the
## `lambda` attribute:
attr(result1$fhat,"lambda")
attr(result2$fhat,"lambda")
