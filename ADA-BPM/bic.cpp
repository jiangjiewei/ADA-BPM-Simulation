#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double bic(arma::mat X, arma::mat M){
int n = X.n_rows;
int p = X.n_cols;
arma::mat F=(arma::eye(p,p)-arma::ones(p,p)/p);
arma::mat FsigmaF=arma::cov(X)*(1-1/n);
arma::mat path=arma::conv_to<arma::mat>::from(M != 0);
double value=n*arma::norm((F*M*FsigmaF + FsigmaF*M * F) / 2 - F, 1) + std::log(n) * arma::accu(arma::trimatu(path, 1));
return value;}