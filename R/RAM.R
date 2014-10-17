RAM <- function(target, N, init.state, init.cov)  {
  d <- length(init.state)
  res=matrix(,N,d)
  res[1,] = rep(0,d)
  t0 = N/2
  alpha_star = 0.234
  S=diag(d)
  
  for (t in 2:N) {
    mean = res[t-1,]  
    U = rnorm(d)
    eta = min(1,d*t^(-2/3))
    x = mean + S %*% U
    alpha = min(1,target(x)/target(res[t-1,]))
    u = runif(1)
    if (u<alpha) {res[t,] = x}
    else {res[t,] = res[t-1,]}
    
    sigma = S %*% (diag(d) + eta*(alpha-alpha_star)* U %*% t(U)/(sum(U^2))) %*% t(S)  
    S = t(chol(sigma))  
  }    
  return(res)
}  