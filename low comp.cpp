#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
void computeMeansK1(const mat& data, int t, double& mean1, double& mean2) {
  // Extract the relevant parts of the data matrix
  vec data_first_part = data.submat(0, 0, 0, t);          // First part of data up to column t
  vec data_second_part = data.submat(0, t + 1, 0, data.n_cols - 1);  // Second part of data from column t+1 to end
  
  // Compute mean1 and mean2
  mean1 = mean(data_first_part);
  mean2 = mean(data_second_part);
}
void computeMeansKGreaterThan1(const mat& data, int t, vec& mean1, vec& mean2) {
  // Extract the relevant parts of the data matrix
  mat data_first_part = data.cols(0, t);          // First part of data up to column t
  mat data_second_part = data.cols(t + 1, data.n_cols - 1);  // Second part of data from column t+1 to end
  
  // Compute means along the rows
  mean1 = mean(data_first_part, 1);  // Mean1 is a vector of means of rows for first part
  mean2 = mean(data_second_part, 1); // Mean2 is a vector of means of rows for second part
}
// [[Rcpp::export]]
void computeStatistics(const mat& data, int n,int k) {
  int sublength = n - 1;
  vec scale = linspace(1, sublength, sublength) / (n * n);
  
  vec substat(sublength);
  
  for (int t = 0; t < sublength; ++t) {
    double mean1, mean2;
    
    if (k == 1) {
      computeMeansK1(data, t, mean1, mean2);
    } else {
      vec mean1_vec, mean2_vec;
      computeMeansKGreaterThan1(data, t, mean1_vec, mean2_vec);
      mean1 = mean(mean1_vec);
      mean2 = mean(mean2_vec);
    }
    
    // Compute inter1 and inter2
    mat inter1(k, t + 1);
    mat inter2(k, n - t - 1);
    
    for (int h = 0; h < k; ++h) {
      inter1.row(h) = cumsum(data.row(h).subvec(0, t)) - (linspace(1, t + 1, t + 1) * mean1);
      inter2.row(h) = cumsum(data.row(h).subvec(t + 1, n - 1)) - (linspace(1, n - t - 1, n - t - 1) * mean2);
    }
    
    // Compute M1 and M2
    mat M1 = (inter1 * trans(inter1)) / (n * n);
    mat M2 = (inter2 * trans(inter2)) / (n * n);
    
    // Compute substat[t]
    // Compute scale * (mean1 - mean2) as a scalar vector
    arma::vec mean_diff = arma::vec(1);  // Create a vector with one element
    mean_diff(0) = scale[t] * (mean1 - mean2);  // Assign the computed value to the first element
    
   // vec mean_diff = scale[t] * (mean1 - mean2);
    substat[t] = n * as_scalar(mean_diff.t() * inv(M1 + M2) * mean_diff);
  }
}
