
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Helper function to extract a row slice from a NumericMatrix
arma::vec extract_row_slice(const arma::mat& data, int row, int start, int end) {
  if (end >= data.n_cols) end = data.n_cols - 1;
  arma::vec slice(end - start + 1);
  for (int i = start; i <= end; i++) {
    slice(i - start) = data(row, i);
  }
  return slice;
}

// Main function to compute W values
// [[Rcpp::export]]
arma::vec crit_val_cpp(int reps, int k, int n, double eps) {
  arma::vec W(reps);
  int sublength = n - 1;
  double n_squared = n * n;
  
  for (int j = 0; j < reps; j++) {
    arma::mat data = arma::randn(k, n);
    arma::vec substat(sublength);
    arma::vec scale(sublength);
    
    // Precompute the scaling factor for each t
    for (int i = 0; i < sublength; i++) {
      scale(i) = (i + 1) * (sublength - i) / n_squared;
    }
    
    for (int t = 0; t < sublength; t++) {
      arma::vec mean1(k), mean2(k);
      
      // Calculate means for each segment
      for (int i = 0; i < k; i++) {
        mean1(i) = mean(extract_row_slice(data, i, 0, t));
        mean2(i) = mean(extract_row_slice(data, i, t + 1, n - 1));
      }
      
      // Calculate inter1 and inter2
      arma::mat inter1(k, t + 1), inter2(k, n - t - 1);
      for (int h = 0; h < k; h++) {
        arma::vec temp1 = cumsum(extract_row_slice(data, h, 0, t)) - (arma::repmat(arma::vec(t + 1) - 1, 1, 1) * mean1(h));
        arma::vec temp2 = cumsum(extract_row_slice(data, h, t + 1, n - 1)) - (arma::repmat(arma::vec(n - t - 1) - 1, 1, 1) * mean2(h));
        inter1.row(h) = temp1.t();
        inter2.row(h) = temp2.t();
      }
      
      // Compute M1 and M2 matrices
      arma::mat M1 = inter1 * inter1.t() / n_squared;
      arma::mat M2 = inter2 * inter2.t() / n_squared;
      
      // Calculate substat[t] using the quadratic form
      arma::vec diff = scale(t) * (mean1 - mean2);
      substat(t) = n * as_scalar(diff.t() * inv(M1 + M2) * diff);
    }
    // Store the maximum of the substatistics adjusted for eps
    int lower_idx = floor(n * eps);
    int upper_idx = ceil(n * (1 - eps)) - 1;
    W(j) = max(abs(substat.subvec(lower_idx, upper_idx)));
  }
  
  // Return W 
  return W;
}
