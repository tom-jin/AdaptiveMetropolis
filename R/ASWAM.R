#' Adaptive Scaling Within Adaptive Metropolis Sampler
#' 
#' @param target The target distribution where the unscaled densities can be 
#' evaluated.
#' @param n The number of samples to produce.
#' @param init.state The initial state of the sampler.
#' @param init.cov The initial covariance of the sampler.
#' @return n samples from the target distribution.
#' @importFrom MASS mvrnorm
#' @export
ASWAM <- function(target, n, init.state, init.cov) {
  d <- length(init.state)
  X <- matrix(NA, nrow = n + 1, ncol = d)
  X[1, ] <- init.state

  lambda <- 1
  sigma <- init.cov
  
  for(i in seq(1, n)) {
    Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = lambda * sigma)
    a <- min(1, target(Y)/target(X[i, ])) 
    
    if(runif(1) < a) 
      X[i+1, ] <- Y
    else
      X[i+1, ] <- X[i, ]
    
    if(d == 1) {
      Xbarold <- mean(X[1:i, ])
      Xbarnew <- mean(X[seq(1, i+1), ])
    } else {
      Xbarold <- colMeans(X[1:i, ])
      Xbarnew <- colMeans(X[seq(1, i+1), ])
    }
    
    lambda <- exp(log(lambda) + a - 0.234)
    sigma <- ((i-1)/i) * sigma + (2.4^2)/(d*i) * ((i * Xbarold %*% t(Xbarold))
                                          - ((i + 1 ) * Xbarnew %*% t(Xbarnew)) 
                                          + (X[i+1, ] %*% t(X[i+1, ])) 
                                          + .Machine$double.eps*diag(d))
  }
  
  return(X[-1, ])
}