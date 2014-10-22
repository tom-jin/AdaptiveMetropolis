validateInput <- function(target, n, init.state, init.cov, burn.in = 0) {
  if(!is.function(target))
    stop("target must be a function that evaluates the target distribution")
  
  if(!is.numeric(n) | length(n) != 1)
    stop("n must be the number of samples to produce")
  
  if(!is.vector(init.state))
    stop("init.state must be a vector containing the initial state")
  
  d <- length(init.state)
  if(d == 1) {
    if(!is.numeric(init.cov) | length(init.cov) != 1)
      stop("init.cov should be a single number")
  } else {
    if(!is.matrix(init.cov))
      stop("init.cov must be a matrix")
    
    if(dim(init.cov)[1] != d | dim(init.cov)[2] != d)
      stop("init.cov has inconsistent dimension")
    
    #TODO: Check for positive definiteness.
    #TODO: Check for symmetry.
  }
  
  if(!is.numeric(burn.in) | burn.in < 0)
    stop("burn.in must be a non-negative number")
}