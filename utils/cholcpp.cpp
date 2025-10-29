#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd cholcpp(Eigen::MatrixXd Sigma_k) {
  Eigen::LLT<Eigen::MatrixXd> llt(Sigma_k);
  Eigen::MatrixXd L = llt.matrixL();
  return L;
}