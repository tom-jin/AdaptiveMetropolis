## ----packages, message=FALSE---------------------------------------------
require(AdaptiveMetropolis) 

## ----AM-code, eval=FALSE-------------------------------------------------
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

## ----AM-trace, cache=TRUE, echo=FALSE, message=FALSE---------------------
plot(X, xlab = "Sample", ylab = "Value")
abline(v = 1000)

## ----AM-variance, cache=TRUE, echo=FALSE, message=FALSE------------------
plot(C, xlab = "Sample", ylab = "Variance")
abline(v = 1000)

## ----ASWAM, cache=TRUE---------------------------------------------------
library(mvtnorm)
data <- ASWAM(function(x) {dmvnorm(x, rep(10, 8), 10*diag(8))}, 
              1000, rep(0, 8), diag(8))
pairs(data)

