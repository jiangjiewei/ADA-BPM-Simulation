#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat Trf(arma::mat M) {
  int p = M.n_rows;  
  double s = arma::accu(M);  
  arma::rowvec v = arma::sum(M, 0);  
  arma::mat Mv = arma::repmat(v, p, 1);  
  arma::mat C = M - Mv / p - Mv.t() / p + std::pow(p, -2) * s * arma::ones(p, p);  
  return C;
}