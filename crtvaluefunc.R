# ===========================================================
# Function to compute critical values for Corollaries 1 and 2
# (based on https://publish.illinois.edu/xshao/files/2013/11/Change-point.txt)
# ===========================================================

# reps:     number of Monte Carlo repetitions
# k:        number of independent Wiener processes to simulate
# n:        number of gridpoints for each Wiener process 
# Quantile: quantiles of interest (= 1 - significance level)
# eps:      \epsilon in the notation of Corollaries 1 and 2
install.packages("Rcpp")
library(Rcpp)
crit.val <- function(reps, k, n, Quantile, eps){
  W         <- rep(0,reps)
  sublength <- n-1
  scale     <- (1:(n-1))/n^2*((n-1):1)
  
  for (j in 1:reps){
    data    <- matrix( rnorm(n*k, 0, 1), k, n )
    substat <- rep(0, sublength)
    for (t in 1:(n-1)){
      if (k==1){
        mean1 <- mean( data[1, 1:t] )
        mean2 <- mean( data[1, (t+1):n] )
      }
      if (k>1){
        mean1 <- apply( matrix(data[,     1:t], ncol = t),   1, mean )
        mean2 <- apply( matrix(data[, (t+1):n], ncol = n-t), 1, mean )
      }
      inter1 <- NULL
      inter2 <- NULL
      for (h in 1:k){
        inter1 <- rbind(inter1, cumsum( data[h, 1:t    ] ) - (1:t)*mean1[h])
        inter2 <- rbind(inter2, cumsum( data[h, n:(t+1)] ) - (1:(n-t))*mean2[h])
      }
      M1 <- inter1 %*% t(inter1)/n^2
      M2 <- inter2 %*% t(inter2)/n^2
      
      substat[t] <- n * matrix(scale[t]*(mean1 - mean2), 1, k) %*% solve(M1+M2) %*% 
                        matrix(scale[t]*(mean1 - mean2), k, 1) 
    }
    W[j] <- max( abs(substat)[ floor(n*eps):ceiling(n*(1-eps)) ] )
    print(j)
  }
  return( quantile( W, Quantile ) )
}

start <- Sys.time()
crit.val(reps=10, k=1, n=10, Quantile = c(0.9,0.95,0.99,0.995,0.999), eps=0.01)
print( Sys.time() - start )


library(RcppArmadillo)
RcppArmadillo::compileArmadilloCode()
install.packages("Rcpp")
library(Rcpp)
Rcpp::sourceCpp('crtvalue.cpp')
initdata <- rnorm(1000, mean = 21, sd = 10)
bootstrap_cpp(initdata)

# Function declaration
bootstrap_r <- function(ds, B = 1000) {
  # Preallocate storage for statistics
  boot_stat <- matrix(NA, nrow = B, ncol = 2)
  # Number of observations
  n <- length(ds)
  # Perform bootstrap
  for(i in seq_len(B)) {
    # Sample initial data
    gen_data <- ds[ sample(n, n, replace=TRUE) ]
    # Calculate sample data mean and SD
    boot_stat[i,] <- c(mean(gen_data),
                       sd(gen_data))
  }
  # Return bootstrap result
  return(boot_stat)
}

set.seed(512)
# Generate data
initdata <- rnorm(1000, mean = 21, sd = 10)
# Set a new _different_ seed for bootstrapping
set.seed(883)
# Perform bootstrap
result_r <- bootstrap_r(initdata)

Rcpp::sourceCpp('try.cpp')
crtvalue(reps=10,k=5,n=10,eps=0.01)
Rcpp::sourceCpp('gggg.cpp')
Rcpp::sourceCpp('Untitled.cpp')
Rcpp::sourceCpp('try.cpp')
Rcpp::sourceCpp('ttt.cpp')

