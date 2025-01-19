#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


//' @title Element-wise Soft Thresholding Function for Matrix
//' @description Soft thresholding function for Matrix
//'
//' @param A Matrix
//' @param a scalar
//' @param diag  Should diagonal elements of the Matrix be thresholding? Default is FALSE.
//' @return Matrix after threholding
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat soft(arma::mat A,double a, int diag=0){
  a=std::abs(a);
  arma::mat C=(A>=a)%(A-a)+(A<=(-a))%(A+a);
if (diag==0){
  arma::mat B=A-arma::diagmat(arma::diagvec(A));
  B=(B>=a)%(B-a)+(B<=(-a))%(B+a);
  B=B+arma::diagmat(arma::diagvec(A));
  C=B;}
  return C;}


//' @title FMF
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

 
// [[Rcpp::export]]
Rcpp::List adabpm(arma::mat X,arma::vec lambda,double err=10^(-5),int maxIter=1000,double gamma=0.0001,int diag=0){
  int n=X.n_rows;
  int p=X.n_cols;
  int nlambda=lambda.size();
  double rho=0.001;
  // double rho=1.01;
  /*Centering int*/
  arma::mat F=(arma::eye(p,p)-arma::ones(p,p)/p);
  arma::mat dX=arma::cov(X)*(1-1/n);
  arma::mat U;
  arma::mat Uv;
  arma::vec eigd0;  
  arma::svd(U, eigd0, Uv, dX.t());
  arma::vec eigd = eigd0.elem(find(eigd0 > 1e-10));
  arma::uvec ind = find(eigd0 > 1e-10);
  arma::mat U1=U.cols(ind)*arma::diagmat(sqrt(eigd/(eigd+2*rho)));
  arma::vec ld=exp(eigd);
  arma::mat D=2*rho/(log(ld*ld.t())+2*rho)+1;
 Rcpp::List Omega_all(nlambda);
 Rcpp::List path(nlambda);
 arma::vec niter(nlambda);
 arma::vec bic(nlambda);
 /*Intialization*/
 arma::mat aZ=arma::eye(p,p);
 arma::mat aU=arma::zeros(p,p);
 arma::mat aX;
 arma::mat L;
 arma::mat L1;
 arma::mat L2;
 arma::mat Z1;
 arma::mat path_mat;
 double lam;
 double ee=1;
 for (int k=0;k<nlambda;++k) {
   lam=lambda(k);
   int i=0;
   while (((i<maxIter)&&(ee>err))||(i==0))
   { Z1=aZ;
     L=F/rho+aZ-aU;
     L=(L+L.t())/2;
     L1=L*U1;
     L2=L1*U1.t();
     L2=F*L2;     
     aX=L-L2-L2.t()+Trf(U1*(D%(U1.t()*L1))*U1.t());
     aZ=soft(aX+aU,lam/rho,diag);
     aU=aU+gamma*(aX-aZ);
     ee=mean(mean(abs(aZ-Z1)));
     i=i+1;
   }
   Omega_all(k)=arma::sp_mat(aZ);
   path(k) = arma::conv_to<arma::mat>::from(aZ != 0);
   niter(k)=i;
   path_mat=Rcpp::as<arma::mat>(path(k));
   bic(k)=n*arma::norm((F*aZ*dX + dX*aZ * F) / 2 - F, 1) + std::log(n) * arma::accu(arma::trimatu(path_mat, 1));
 }
 return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                           Rcpp::Named("lambda") =lambda,
                           Rcpp::Named("niter") =niter,
                           Rcpp::Named("path") =path,
                           Rcpp::Named("bic") =bic); }
