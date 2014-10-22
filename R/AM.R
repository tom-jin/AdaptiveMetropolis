#' Adaptive Metropolis Sampler
#' 
#' @param target The target distribution where the unscaled densities can be 
#' evaluated.
#' @param n The number of samples to produce.
#' @param init.state The initial state of the sampler.
#' @param init.cov The initial covariance of the sampler.
#' @return n samples from the target distribution.
#' @importFrom MASS mvrnorm
#' @export
#' @examples
#' data <- AM(function(x) {dunif(x, 0,100)}, 1000, 0, 1)
#' plot(data)
#' data <- AM(function(x) {dunif(x, c(0,0), c(100,100))}, 1000, c(0,0), diag(2))
#' plot(data)
AM <- function(target, n, init.state, init.cov) {
  d <- length(init.state)
  X <- matrix(NA, nrow = n + 1, ncol = d)
  X[1, ] <- init.state
  C <- init.cov
  
  for(i in seq(1, n)) {
    # Propose a move
    Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = C)
    if(is.nan(Y[1])) stop ("Balls.")
    # Calculate acceptance ratio
    a <- min(1, target(Y)/target(X[i, ])) 
    
    if(runif(1) < a) 
      X[i+1, ] <- Y # Accept
    else
      X[i+1, ] <- X[i, ] # Reject
    
    Xbarold <- .colMeans(X[1:i, ], i, d)
    Xbarnew <- .colMeans(X[seq(1, i+1), ], i+1, d)
    
    # Update covariance
    C <- ((i-1)/i) * C + (2.4^2)/(d*i) * ((i * Xbarold %*% t(Xbarold))
                                          - ((i + 1 ) * Xbarnew %*% t(Xbarnew)) 
                                          + (X[i+1, ] %*% t(X[i+1, ])) 
                                          + .Machine$double.eps*100*diag(d))
    #C <- cov(X[seq(1, i+1), ])
  }
  return(X[-1, ])
}