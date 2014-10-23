#' Adaptive Metropolis-Within-Gibbs
#' 
#' Adaptive Metropolis-Hastings algorithm within the Gibbs algorithm, applied to
#' a specific hierachical model.
#' 
#' Update parameter one by one using a normal proposal, accept it with 
#' probability depending on the simulated path.
#' 
#' @param adapt whether to use adapted MH or normal MH algorithm
#' @param N number of iterations
#' @param K number of latent variables
#' @return A list containing simulated parameters paths (theta) and 
#' acceptance/rejection paths (acceptance).
#' @references Roberts, G. O., Rosenthal, J. S. (2006), 
#' \emph{Examples of Adaptive MCMC}
#' @author Xiaoyu Lu
#' @export
#' @examples
#' out_noadapt <- AMWG(adapt = 0, 1000, 100)
#' theta_noadapt <- out_noadapt[[1]]
#' out_adapt <- AMWG(adapt = 1, 1000, 100)
#' theta_adapt <- out_adapt[[1]]
#' plot(theta_noadapt[,1])
#' points(theta_adapt[,1],col=2)
AMWG <- function(adapt, N, K) {
  acceptance <- matrix(0, N, K+3)
  
  # Initial condition
  theta <- matrix(, N, K+3)
  theta[1, ] <- rep(100, K+3)
  alpha <- c()
  ls <- rep(0, K+3)
  
  #Simulated data
  r <- sample(5:500, size = K, replace = TRUE)
  Y <- list()
  for (i in 1:K) {
    Y[[i]] <- rnorm(r[i]) * 10 + (i-1) 
  }
  
  #Gibbs sampler
  for(n in 1:(N-1))  {
    delta <- min(0.001, n^(-1/2))
    
    if(n > 50) {
      # adjust the proposal deviation according to the proportion of accepted proposals up to time n
      index <- which(colSums(acceptance)/n > 0.44)
      ls[index] <- ls[index]+delta*adapt
      ls[setdiff(1:(K+3), index)] <- ls[setdiff(1:(K+3), index)]-delta*adapt
    }
    
    #bound the standard deivation
    for(k in 1:(K+3)) {
      ls[k] <- max(min(ls[k], 3), -3)
    }
    
    #propose new moves
    theta_new <- rnorm(K+3) * exp(ls) + theta[n, ]
    u <- runif(K+3)
    # append A, V, mu to theta and calculate acceptance probabilities
    for (i in 1:K) {    
      alpha[i] <- min(1, (1+((theta[n, i]-theta[n, K+3])/theta[n, K+1])^2)/
                        (1+((theta_new[i]-theta[n, K+3])/theta[n, K+1])^2) *
                        exp(1/(2*theta[n, K+2])*sum((Y[[i]]-theta[n, i])^2-
                                                      (Y[[i]]-theta_new[i])^2)))       
    }
    
    alpha[K+1] <- min(1, exp(-1/theta_new[K+1]+1/theta[n, K+1])*
                        (theta[n, K+1]/theta_new[K+1])^2 * 
                        ifelse(theta_new[K+1] >= 0, 1, 0) * 
                        prod(1+((theta[n, 1:K]-theta[n, K+3])/theta[n, K+1])^2)/
                        prod(1+((theta_new[1:K]-theta[n, K+3])/theta[n, K+1])^2), 
                      na.rm = T)
    
    alpha[K+2] <- min(1, (theta[n, K+2]/theta_new[K+2])^(sum(r)/2+2) * 
                        ifelse(theta_new[K+2]>=0, 1, 0) *
                        exp((1/theta[n, K+2]-1/theta_new[K+2])*
                              (1+sum((unlist(Y)-rep(theta[n, 1:K], r))^2)/2)), 
                      na.rm=T)
    
    alpha[K+3] <- min(1, exp(-theta_new[K+3]^2/2+theta[n, K+3]^2/2) 
                      * prod(1+((theta[n, 1:K]-theta[n, K+3])/theta[n, K+1])^2)/
                        prod(1+((theta_new[1:K]-theta[n, K+3])/theta[n, K+1])^2), 
                      na.rm = T)
    
    index1 <- which(u < alpha)
    index2 <- setdiff(1:(K+3), index1)
    theta[n+1, index1] <- theta_new[index1]
    theta[n+1, index2] <- theta[n, index2]
    acceptance[n, which(u < alpha)] <- 1
    
  }
  return(list(theta, acceptance))
}