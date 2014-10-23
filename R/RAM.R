#' Robust Adaptive Metropolis Sampler
#' 
#' RAM function takes inputs including target distribution, number of 
#' iterations, initial state and initial proposal covariance matrix
#' 
#' @param target The target distribution where the unscaled densities can be 
#' evaluated.
#' @param n The number of samples to produce.
#' @param init.state The initial state of the sampler.
#' @param init.cov The initial covariance of the sampler.
#' @return n samples from the target distribution.
#' @export
RAM <- function(target, N, init.state, init.cov)  {
  
  d <- length(init.state)  # dimension of the space of interest
  res <- matrix(, N, d)         # set up the result matrix, where each row n represents the parameter of length d at time n 
  res[1, ] <- init.state     # store the initial state as the first result
  t0 <- N/2                 # burn-in period is set to be N/2, where N is the number of iterations
  alpha_star <- 0.234       # target mean acceptance rate is set to be 0.234 according to some popular theoretical result
  S <- diag(d)                # initial covariance is set to be the identity matrix, i.e. components are independent
  
  for (t in 2:N) {
    eta <- min(1, d*t^(-2/3))     # set eta
    
    mean <- res[t-1, ]       # set the mean for the proposal to be the current state
    U <- rnorm(d)           # simulate a d-dimensional MVN variable
    x <- mean + S %*% U          # calculate the proposed point
    alpha <- min(1, target(x)/target(res[t-1, ]))     #calculate acceptance probability
    u <- runif(1)                # simulate a uniform random variable and accept/reject the proposal accordingly.
    if (u<alpha) 
      res[t, ] <- x
    else 
      res[t, ] <- res[t-1, ]
    
    sigma <- S %*% (diag(d) + eta*(alpha-alpha_star)* U %*% t(U)/(sum(U^2))) %*% t(S)  #update the covariance matrix
    S <- t(chol(sigma))  
  }    
  return(res)         # output the simulated path
}  


