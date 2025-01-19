BIC<-function(X,theta){
  n<-nrow(X)
  p<-ncol(X)
  F <- diag(p) - matrix(1,p,p)/p
  FsigmaF<-cov(X)*(1-1/n)
  path <- 1*(theta!=0)
  bic <- n*norm((F%*%theta%*%FsigmaF+FsigmaF%*%theta%*%F)/2-F,"1")+
    log(n)*sum(path[upper.tri(path)])
  return(bic)
}
