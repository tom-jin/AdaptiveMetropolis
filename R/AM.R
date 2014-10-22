#' Adaptive Metropolis Sampler
#' 
#' An implementation of the Adaptive Metropolis algorithm as proposed by Haario 
#' et al (2001). 
#' The covariance matrix of the proposal distribution is adapted with current 
#' knowledge of the target distrubution.
#' 
#' The proposal distribution is a multivariate normal distribution of the 
#' appropriate dimension. The covariance of this distribution is computed from
#' previous samples.
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
#' @param target A function that evaluates the unscaled densities of the target 
#' distribution. 
#' @param n The number of samples to produce.
#' @param init.state The initial state of the sampler.
#' @param init.cov The covariance of the sampler during burn in.
#' @param burn.in (Optional) Number of iterations to burn in for. Default: 1000
#' @return \code{n} samples from the target distribution.
#' @references Haario, H., Saksman, E., Tamminen, J. (2001) 
#' \emph{An Adaptive Metropolis Algorithm.} Bernoulli, Vol. 7, No. 2.
#' @author Tom Jin <sjin@@stats.ox.ac.uk>
#' @importFrom MASS mvrnorm
#' @export
#' @examples
#' # 1D Normal Distribution
#' data <- AM(function(x) {dnorm(x, 0,100)}, 5000, 0, 1)
#' plot(data)
#' 
#' # 2D Uniform Distribution
#' data <- AM(function(x) {dunif(x, c(0,0), c(100,100))}, 5000, c(0,0), diag(2))
#' plot(data)
AM <- function(target, n, init.state, init.cov, burn.in = 1000) {
  # Input validation
  validateInput(target, n, init.state, init.cov, burn.in)
  
  # Init variables
  d <- length(init.state)
  X <- matrix(NA, nrow = n + burn.in + 1, ncol = d)
  X[1, ] <- init.state
  
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
    C <- var(X[seq(1, burn.in), ])
  } else {
    C <- init.cov
  }
  
  # Adapted Metropolis phase
  for(i in seq(1, n) + burn.in) {
    # Propose a move
    Y <- mvrnorm(n = 1, mu = X[i, ], Sigma = C)

    # Calculate acceptance ratio
    a <- min(1, target(Y)/target(X[i, ])) 
    
    if(runif(1) < a) 
      X[i+1, ] <- Y # Accept
    else
      X[i+1, ] <- X[i, ] # Reject
    
    Xbarold <- .colMeans(X[1:i, ], i, d)
    Xbarnew <- .colMeans(X[seq(1, i+1), ], i+1, d)
    
    # Update covariance
    C <- ((i-1)/i) * C + 
      (2.4^2)/(d*i) * ((i * Xbarold %*% t(Xbarold))
                       - ((i + 1) * Xbarnew %*% t(Xbarnew)) 
                       + (X[i+1, ] %*% t(X[i+1, ])) 
                       + .Machine$double.eps*diag(d))
  }
  return(X[-seq(1, burn.in + 1), ])
}
