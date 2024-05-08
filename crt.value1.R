library(RCPP)
Rcpp::sourceCpp('arma.cpp')

crit.val2 <- function(reps, k, n, Quantile, eps){
  vcrt <- crit_val_cpp(reps, k, n, eps)
  return( quantile( vcrt, Quantile ) )
  }
start <- Sys.time()
crit.val2(reps=10, k=1, n=10, Quantile = c(0.9,0.95,0.99,0.995,0.999), eps=0.01)
print( Sys.time() - start )





