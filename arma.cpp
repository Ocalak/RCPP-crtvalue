#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

// Helper function to extract a row slice from a NumericMatrix
Rcpp::NumericVector extract_row_slice(const Rcpp::NumericMatrix& data, int row, int start, int end) {
  if (end >= data.ncol()) end = data.ncol() - 1;
  Rcpp::NumericVector slice(end - start + 1);
  for (int i = start; i <= end; i++) {
    slice[i - start] = data(row, i);
  }
  return slice;
}

// Main function to compute critical values
// [[Rcpp::export]]
Rcpp::NumericVector getW(int reps, int k, int n, double eps) {
  Rcpp::NumericVector W(reps);
  int sublength = n - 1;
  double n_squared = n * n;
  
  for (int j = 0; j < reps; j++) {
    Rcpp::NumericMatrix data(k, n, Rcpp::rnorm(n * k).begin());
    Rcpp::NumericVector substat(sublength);
    Rcpp::NumericVector scale(sublength);
    
    // Precompute the scaling factor for each t
    for (int i = 0; i < sublength; i++) {
      scale[i] = (i + 1) * (sublength - i) / n_squared;
    }
    
    for (int t = 0; t < sublength; t++) {
      Rcpp::NumericVector mean1(k), mean2(k);
      
      // Calculate means for each segment
      for (int i = 0; i < k; i++) {
        mean1[i] = mean(extract_row_slice(data, i, 0, t));
        mean2[i] = mean(extract_row_slice(data, i, t + 1, n - 1));
      }
      
      // Calculate inter1 and inter2
      Rcpp::NumericMatrix inter1(k, t + 1), inter2(k, n - t - 1);
      for (int h = 0; h < k; h++) {
        Rcpp::NumericVector temp1 = Rcpp::as<Rcpp::NumericVector>(wrap(cumsum(extract_row_slice(data, h, 0, t)))) - Rcpp::as<Rcpp::NumericVector>(wrap(Rcpp::seq_len(t + 1) - 1)) * mean1[h];
        
        Rcpp::NumericVector temp2 = Rcpp::as<Rcpp::NumericVector>(wrap(cumsum(extract_row_slice(data, h, t + 1, n - 1)))) - Rcpp::as<Rcpp::NumericVector>(wrap(Rcpp::seq_len(n - t - 1) - 1)) * mean2[h];
        
        std::copy(temp1.begin(), temp1.end(), inter1(h,Rcpp::_).begin());
        std::copy(temp2.begin(), temp2.end(), inter2(h,Rcpp::_).begin());
      }
      
      // Compute M1 and M2 matrices
      Rcpp::NumericMatrix M1 = inter1 * inter1.t() / n_squared;// Issue here. Error1: t() function is not exist. I use arma:trans it didnt work as well. 
      Rcpp::NumericMatrix M2 = inter2 * inter2.t() / n_squared;// Same issue as well. 
      
      // Calculate substat[t] using the quadratic form
      Rcpp::NumericVector diff = scale[t] * (mean1 - mean2);
      substat[t] = n * crossprod(diff,solve(M1 + M2, diff)); // Last issue. Cannot use solve() function.  arma::solve doesnt work as well. 
    }
    
    // Store the maximum of the substatistics adjusted for eps
    int lower_idx = floor(n * eps);
    int upper_idx = ceil(n * (1 - eps)) - 1;
    W[j] = max(abs(substat[Rcpp::Range(lower_idx, upper_idx)]));
  }
  
  // Return W 
  return W;
}



