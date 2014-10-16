#' Adaptive Metropolis Sampler
#' 
#' @param target The target distribution where the unscaled densities can be 
#' evaluated.
#' @param n The number of samples to produce.
#' @param init.state The initial state of the sampler.
#' @param init.cov The initial covariance of the sampler.
#' @param init.period The number of samples to use init.cov instead of adapting 
#' the covariance.
#' @return n samples from the target distribution.
#' @examples
#' data <- AM(function(x) {dunif(x, 0,100)}, 1000, 0, 1, 10)
#' plot(data)
#' data <- AM(function(x) {dunif(x, c(0,0), c(100,100))}, 1000, c(0,0), diag(2), 10)
#' plot(data)
AM <- function(target, n, init.state, init.cov, init.period = 1000) {
  d <- length(init.state)
  X <- matrix(NA, nrow = n + 1, ncol = d)
  X[1, ] <- init.state
  
  for(i in 1:init.period) {
    Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = init.cov)
    a <- min(1, target(Y)/target(X[i, ])) 
    
    if(runif(1) < a) 
      X[i+1, ] <- Y
    else
      X[i+1, ] <- X[i, ]
  }
  C <- var(X[1:init.period, ])
  for(i in seq(init.period+1, n)) {
    Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = C)
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
    C <- ((i-1)/i) * C + (2.4^2)/(d*i) * ((i * Xbarold %*% t(Xbarold))
                                          - ((i + 1 ) * Xbarnew %*% t(Xbarnew)) 
                                          + (X[i+1, ] %*% t(X[i+1, ])) 
                                          + .Machine$double.eps*diag(d))
  }
  
  return(X[-1, ])
}