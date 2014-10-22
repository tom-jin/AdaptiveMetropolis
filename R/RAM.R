#' Robust Adaptive Metropolis Sampler
#' 
#' @param target The target distribution where the unscaled densities can be 
#' evaluated.
#' @param n The number of samples to produce.
#' @param init.state The initial state of the sampler.
#' @param init.cov The initial covariance of the sampler.
#' @return n samples from the target distribution.
#' @export
RAM <- function(target, n, init.state, init.cov)  {
  # Input validation
  validateInput(target, n, init.state, init.cov)
  
  # Init variables
  d <- length(init.state)
  res <- matrix(,n,d)
  res[1,] <- rep(0,d)
  t0 <- n/2
  alpha_star <- 0.234
  S <- diag(d)
  
  for (t in 2:n) {
    mean <- res[t-1,]  
    U <- rnorm(d)
    eta <- min(1,d*t^(-2/3))
    x <- as.vector(mean + S %*% U)
    alpha <- min(1,target(x)/target(res[t-1,]))
    u <- runif(1)
    
    if (u<alpha) 
      res[t,] <- x
    else 
      res[t,] <- res[t-1,]
    
    sigma <- S %*% (diag(d) + eta*(alpha-alpha_star)* U %*% t(U)/(sum(U^2))) %*% t(S)  
    S <- t(chol(sigma))  
  }    
  return(res)
}  