#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd solve_linear_equation(Eigen::MatrixXd Sigma, Eigen::MatrixXd diag_s) {
  Eigen::MatrixXd x = Sigma.ldlt().solve(diag_s);
  return x;
}