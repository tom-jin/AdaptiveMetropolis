## ----packages, message=FALSE, cache=TRUE---------------------------------
require(AdaptiveMetropolis) 

## ----AM-code, eval=FALSE, cache=TRUE-------------------------------------
#  data <- AM(function(x) {dnorm(x, 100,100)}, 5000, 0, 1)
#  plot(data)

## ----AM-plots, cache=TRUE, echo=FALSE, message=FALSE---------------------
require(MASS)
target <- function(x) {dnorm(x, 100,100)}
n <- 5000
init.state <- 0
init.cov <- 1
burn.in <- 1000

# Init variables
  d <- length(init.state)
  X <- matrix(NA, nrow = n + burn.in + 1, ncol = d)
  X[1, ] <- init.state
  C <- rep(init.cov, n + burn.in + 1)

  # Vanilla Metropolis phase. (Burn in)
  if (burn.in > 0)
    for(i in seq(1, burn.in)) {
      # Propose a move
      Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = init.cov)
      
      # Calculate acceptance ratio
      a <- min(1, target(Y)/target(X[i, ])) 
      
      if(runif(1) < a) 
        X[i+1, ] <- Y # Accept
      else
        X[i+1, ] <- X[i, ] # Reject
    }
  
  # Initialise covariance
  if(burn.in > 1) {
    C[burn.in+1] <- var(X[seq(1, burn.in), ])
  }

  # Adapted Metropolis phase
  for(i in seq(1, n) + burn.in) {
    # Propose a move
    Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = C[i])

    # Calculate acceptance ratio
    a <- min(1, target(Y)/target(X[i, ])) 
    
    if(runif(1) < a) 
      X[i+1, ] <- Y # Accept
    else
      X[i+1, ] <- X[i, ] # Reject
    
    Xbarold <- .colMeans(X[1:i, ], i, d)
    Xbarnew <- .colMeans(X[seq(1, i+1), ], i+1, d)
    
    # Update covariance
    C[i+1] <- ((i-1)/i) * C[i] + 
      (2.4^2)/(d*i) * ((i * Xbarold %*% t(Xbarold))
                       - ((i + 1) * Xbarnew %*% t(Xbarnew)) 
                       + (X[i+1, ] %*% t(X[i+1, ])) 
                       + .Machine$double.eps*diag(d))
  }
#par(mfrow=c(1,2))

## ----AM-trace, cache=TRUE, echo=FALSE, message=FALSE, fig.height=4-------
plot(X, xlab = "Sample", ylab = "Value")
abline(v = 1000)

## ----AM-variance, cache=TRUE, echo=FALSE, message=FALSE, fig.height=4----
plot(C, xlab = "Sample", ylab = "Variance")
abline(v = 1000)

## ----ASWAM, cache=TRUE---------------------------------------------------
library(mvtnorm)
data <- ASWAM(function(x) {dmvnorm(x, rep(10, 8), 10*diag(8))}, 
              1000, rep(0, 8), diag(8))
pairs(data)

## ----chunk1, cache=TRUE--------------------------------------------------
library(mvtnorm)  
res <- RAM(function(x) {dmvnorm(x, rep(0, 3), diag(3))}, 10000, 
           c(0,0,0), diag(3))

## ----chunk2, cache=TRUE, echo=FALSE, fig.height=3------------------------
plot(res[,1],type = "l",xlab="Iteration",main="Traceplot")

## ----chunk2a, cache=TRUE, echo=FALSE, fig.height=3-----------------------
hist(res[,1],main="Histogram")

## ----chunk3, cache=TRUE--------------------------------------------------
out_unadapt <- AMWG(adapt = 0, 1000, 500)
theta_unadapt <- out_unadapt[[1]]
acceptance_unadapt <- out_unadapt[[2]]
out_adapt <- AMWG(adapt = 1, 1000, 500)
theta_adapt <- out_adapt[[1]]
acceptance_adapt <- out_adapt[[2]]

## ----chunk, cache=TRUE, echo=FALSE, fig.height=5-------------------------
plot(theta_unadapt[,1],type="l",xlab="Batch",ylab=expression(theta[1]),main=expression('Traceplots of '*theta[1]*' using Metropolis within Gibbs'))
lines(theta_adapt[,1],col=2)
legend("topright",c("unadapted","adapted"),col=c(1,2),lwd=c(2,2))

## ----chunk4, cache=TRUE, echo=FALSE, fig.height=5------------------------
par(mfrow=c(2,1))
hist(colMeans(acceptance_unadapt),main="Histogram of acceptance ratio for K+3 parameters (unadapted)",xlab="",breaks=50)
hist(colMeans(acceptance_adapt),main="Histogram of acceptance ratio for K+3 parameters (adapted)",xlab="",breaks=50)

## ----echo=FALSE, cache=TRUE, fig.height=3--------------------------------
library(mnormt)
N=10000
init.state = rnorm(2)    #random initialization
d = length(init.state)
init.cov = diag(d)
mean = c(1,2)
sigma = matrix(c(0.2,0.1,0.1,0.8),2,2)
target = function(x) {
  return(dmt(t(x),mean,sigma,1))
}

resAM <- AM(target, N, init.state, init.cov)
resASWAM <- ASWAM(target, N, init.state, init.cov)
resRAM <- RAM(target, N, init.state, init.cov)

par(mfrow=c(1,3))
par(pty="s")
plot(resAM[,1],resAM[,2])
plot(resASWAM[,1],resASWAM[,2])
plot(resRAM[,1],resRAM[,2])

