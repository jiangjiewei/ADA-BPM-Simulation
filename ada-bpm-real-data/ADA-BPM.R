#' Efficient admm algorithm via the QUAdratic Loss (EQUAL) for precision matrix estimation
#' @param X data matrix of dimension n*p.
#' @param type Should the loss function be symmetric? Default is TRUE.
#' @param sdiag Should diagonal of inverse covariance be penalized? Default is FALSE.
#' @param lambda user supplied tuning parameter; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. 
#' Default is sqrt(log(p)/n).
#' @param err the precision used to stop the convergence. Default is 1e-5. 
#' Iterations stop when average absolute parameter change is less than \code{err}.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A list with components
#' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
#' \item{lambda}{the used lambda for the solution path.}
#' \item{niter}{the number of iterations for each element of lambda.}
ADA_BPM<-function(X,sdiag=FALSE,lambda=NULL,lambda.min.ratio=0.001, nlambda=50,err=10^(-5),maxIter =1000,gamma=0.0001)
{p=ncol(X);
n=nrow(X);
FsigmaF<-cov(X)*(1-1/n);
F=diag(p)-matrix(1,p,p)/p;
if (is.null(lambda)){
  lam.max <- max(max(FsigmaF - diag(p)), -min(FsigmaF - diag(p)));
  lam.min <- lambda.min.ratio * lam.max;
  lambda <- exp(seq(log(lam.max), log(lam.min), length = nlambda));
}

obj<-adabpm(X,lambda =lambda,diag=as.numeric(sdiag),err=err,maxIter =maxIter,gamma=gamma)
return(obj)
}  