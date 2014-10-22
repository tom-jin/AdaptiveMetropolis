#' Adaptive Scaling Within Adaptive Metropolis Sampler
#' 
#' An implementation of the Adaptive Scaling within Adaptive Metropolis algorithm
#' as published by Andrieu and Thoms (2008).
#' The adapted covariance matrix which encodes knowledge of the second moments 
#' of the target distribution is scaled to achieve a target acceptance ratio.
#' 
#' The proposal distribution is a multivariate normal distribution of the 
#' appropriate dimension. The covariance of this distribution is computed from
#' previous samples and then scaled to aim for an acceptance rate of 0.234. This 
#' is done to stike a balence between a high acceptance rate with poor mixing 
#' and a low acceptance rate with periods of stationarity. Using diffusion 
#' limits an acceptance rate of 0.234 has been found to be optimal.
#' 
#' This implementation does not have a burn in period and it is expected that 
#' the user will select an appropriate number of samples to discard. Samples 
#' generated whilst little is know about the target distribution may be highly 
#' correlated and localised.
#' 
#' The dimension is detected using the length of the \code{init.state}
#' parameter.
#' 
#' Adaptive samplers do not attempt to solve the problem of multi modal 
#' distributions. If the target distribution is known to have multiple modes, 
#' particually if they are far apart, more advanced techniques are recomended 
#' such as nesting adaptive Metropolis within parallel tempering.
#' 
#' @param target The target distribution where the unscaled densities can be 
#' evaluated.
#' @param n The number of samples to produce.
#' @param init.state The initial state of the sampler.
#' @param init.cov The covariance of the sampler during burn in.
#' @param burn.in (Optional) Number of iterations to burn in for. Default: 1000
#' @return \code{n} samples from the target distribution.
#' @references Andrieu, C., Thoms, J. (2008) \emph{A tutorial on adaptive MCMC.} 
#' Statistics and Computing, Vol. 18, Issue 4
#' @author Tom Jin <sjin@@stats.ox.ac.uk>
#' @importFrom MASS mvrnorm
#' @export
#' @examples 
#' # 4D Multivariate normal
#' library(mvtnorm)
#' data <- ASWAM(function(x) {dmvnorm(x, rep(10, 4), 10*diag(4))}, 8000, rep(0, 4), diag(4))
#' pairs(data)
ASWAM <- function(target, n, init.state, init.cov, burn.in = 1000) {
  # Input validation
  validateInput(target, n, init.state, init.cov, burn.in)
  
  # Init variables
  d <- length(init.state)
  X <- matrix(NA, nrow = n + burn.in + 1, ncol = d)
  X[1, ] <- init.state
  lambda <- 1
  
  
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
    sigma <- var(X[seq(1, burn.in), ])
  } else {
    sigma <- init.cov
  }
  
  for(i in seq(1, n) + burn.in) {
    # Propose a move
    Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = lambda * sigma)
    
    # Calculate acceptance ratio
    a <- min(1, target(Y)/target(X[i, ])) 
    
    if(runif(1) < a) 
      X[i+1, ] <- Y # Accept
    else
      X[i+1, ] <- X[i, ] # Reject
    
    Xbarold <- .colMeans(X[1:i, ], i, d)
    Xbarnew <- .colMeans(X[seq(1, i+1), ], i+1, d)
    
    # Update scaling ratio
    lambda <- exp(log(lambda) + a - 0.234)
    
    # Update covariance
    sigma <- ((i-1)/i) * sigma + (2.4^2)/(d*i) * ((i * Xbarold %*% t(Xbarold))
                                          - ((i + 1 ) * Xbarnew %*% t(Xbarnew)) 
                                          + (X[i+1, ] %*% t(X[i+1, ])) 
                                          + .Machine$double.eps*diag(d))
  }
  return(X[-seq(1, burn.in + 1), ])
}
