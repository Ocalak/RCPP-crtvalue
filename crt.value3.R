library(Rcpp)
Rcpp::sourceCpp('helper.cpp')
crit.val3 <- function(reps, k, n, Quantile, eps){
  W         <- rep(0,reps)
  sublength <- n-1
  scale     <- (1:(n-1))/n^2*((n-1):1)
  n_squared <- n^2
  for (j in 1:reps){
    data    <- matrix( rnorm(n*k, 0, 1), k, n )
    substat <- rep(0, sublength)
    
    xx <- csbst(reps,sublength,k,n,eps,n_squared,data,
    scale,substat)
    
    W[j] <- max( abs(xx)[ floor(n*eps):ceiling(n*(1-eps)) ] )
    print(j)
    print(Sys.time())
  }
  return( quantile( W, Quantile ) )
}
start <- Sys.time()
crit.val3(reps=100, k=1, n=1000, Quantile = c(0.9,0.95,0.99,0.995,0.999), eps=0.01)
print( Sys.time() - start )
