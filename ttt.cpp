#include <Rcpp.h>
using namespace Rcpp;

// Helper function to extract a row slice from a NumericMatrix
NumericVector extract_row_slice(const NumericMatrix& data, int row, int start, int end) {
  if (end >= data.ncol()) end = data.ncol() - 1;
  NumericVector slice(end - start + 1);
  for (int i = start; i <= end; i++) {
    slice[i - start] = data(row, i);
  }
  return slice;
}

// Main function to compute critical values
// [[Rcpp::export]]
NumericVector crit_val_cpp(int reps, int k, int n, double Quantile, double eps) {
  NumericVector W(reps);
  int sublength = n - 1;
  double n_squared = n * n;
  
  for (int j = 0; j < reps; j++) {
    NumericMatrix data(k, n, Rcpp::rnorm(n * k).begin());
    NumericVector substat(sublength);
    NumericVector scale(sublength);
    
    // Precompute the scaling factor for each t
    for (int i = 0; i < sublength; i++) {
      scale[i] = (i + 1) * (sublength - i) / n_squared;
    }
    
    for (int t = 0; t < sublength; t++) {
      NumericVector mean1(k), mean2(k);
      
      // Calculate means for each segment
      for (int i = 0; i < k; i++) {
        mean1[i] = mean(extract_row_slice(data, i, 0, t));
        mean2[i] = mean(extract_row_slice(data, i, t + 1, n - 1));
      }
      
      // Calculate inter1 and inter2
      NumericMatrix inter1(k, t + 1), inter2(k, n - t - 1);
      for (int h = 0; h < k; h++) {
        NumericVector temp1 = cumsum(extract_row_slice(data, h, 0, t)) - (seq_len(t + 1) - 1) * mean1[h];
        NumericVector temp2 = cumsum(extract_row_slice(data, h, t + 1, n - 1)) - (seq_len(n - t - 1) - 1) * mean2[h];
        std::copy(temp1.begin(), temp1.end(), inter1(h, _).begin());
        std::copy(temp2.begin(), temp2.end(), inter2(h, _).begin());
      }
      
      // Compute M1 and M2 matrices
      NumericMatrix M1 = inter1.t() * inter1 / n_squared;
      NumericMatrix M2 = inter2.t() * inter2 / n_squared;
      
      // Calculate substat[t] using the quadratic form
      NumericVector diff = scale[t] * (mean1 - mean2);
      substat[t] = n * crossprod(diff, solve(M1 + M2, diff));
    }
    
    // Store the maximum of the substatistics adjusted for eps
    int lower_idx = floor(n * eps);
    int upper_idx = ceil(n * (1 - eps)) - 1;
    W[j] = max(abs(substat[Range(lower_idx, upper_idx)]));
  }
  
  // Return quantile of computed statistics
  return W;
}
