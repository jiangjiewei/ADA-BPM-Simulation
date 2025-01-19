#------------------------------------- Set path ---------------------------------------#
rm(list=ls())
setwd('~/ada-bpm')

#----------------------------------- Load packages ------------------------------------#
library(ggplot2)
library(huge)
library(mvtnorm)
library(PRIMAL)
library(RColorBrewer)
library('MASS')
library('Rcpp')
library('Matrix')

#---------------------------------- Import functions ----------------------------------#
source("Care.R")
source("CD-trace.R")
source("gCoda.R")
source("SPIEC-EASI.R")
source("Summary_roc.R")
sourceCpp("Trf.cpp")
sourceCpp('adabpm.cpp')
source("ada-bpm.R")


ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}

#------------------------- Simulations for compositional data -------------------------#
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#

#------------------------------ (Band graph, model (a)) -------------------------------#
set.seed(12345)
n = 200
p = 50  # p = 50，100, 200, 400,800
band_graph = huge.generator(n, p, graph = "band", g = 2, verbose = FALSE)
band_am = as.matrix(band_graph$theta)
band_graph1 = huge.generator(n, p, graph = "band", g = 1, verbose = FALSE)
band_am1 = as.matrix(band_graph1$theta)
Omega_0 = band_am * 0.5 + band_am1 * 0.3 + diag(runif(p, 1, 2))
Omega_0 = (abs(min(eigen(Omega_0)$values)) + 0.01) * diag(p) + Omega_0
Sigma_0 = solve(Omega_0)
ratio = c(1/2, 7/12, 2/3, 3/4) 
Care_result = NULL; ada_bpm_result = NULL; CD_result = NULL
gCoda_result = NULL; Spiec_easi_result = NULL



for(i in 1 : 100){
  Y = rmvnorm(n, rep(0, p), Sigma_0)
  W = exp(Y)
  X = W / rowSums(W)
  G = diag(p) - matrix(1, p, p) / p
  Z = log(X) %*% G
  Nlambda = 50
  FsigmaF<-cov(Z)*(1-1/n);
  lam.max <- max(max(FsigmaF - diag(p)), -min(FsigmaF - diag(p)));
  lambda.min.ratio4=0.001;
  lam.min4 <- lambda.min.ratio4 * lam.max;
  lam4 <- exp(seq(log(lam.max), log(lam.min4), length = Nlambda));
  

  #----- Care -----#
  Care_tt<-ttime(Care_output <- Care_est(X, nlambda = Nlambda, lambda_min = 0.01, ratio = ratio[p/50], nsplit = 5))
  Care_Omega_hat = Care_output$Omega_hat
  Care_result = cbind(c(norm(Care_Omega_hat - Omega_0, "2"),
                        norm(Care_Omega_hat - Omega_0, "1"),
                        norm(Care_Omega_hat - Omega_0, "F"),
                        ROC_PR(Care_Omega_hat, Omega_0),Care_tt,min(eigen(Care_Omega_hat)$values)), Care_result)
  
  
  #----- ada_bpm -----#
  ada_bpm_tt<-ttime(ada_bpm_output <- ADA_BPM(Z,lambda =lam4, gamma=0.98))
  ada_bpm_Omega_hat=ada_bpm_output$Omega[which.min(ada_bpm_output$bic)][[1]]
  ada_bpm_Omega_hat=as.matrix(ada_bpm_Omega_hat)
  ada_bpm_result = cbind(c(Matrix::norm(ada_bpm_Omega_hat - Omega_0, "2"),
                          norm(ada_bpm_Omega_hat - Omega_0, "1"),
                          norm(ada_bpm_Omega_hat - Omega_0, "F"),
                          ROC_PR(ada_bpm_Omega_hat, Omega_0), ada_bpm_tt,min(eigen(ada_bpm_Omega_hat)$values)), ada_bpm_result)
  
  
  #----- CD-trace -----#
  lambda = seq(2 * (log(p) / n) ^ (1 / 2), 10 * (log(p) / n) ^ (1 / 2), length.out = Nlambda)
  CD_tt<-ttime(CD_output <- CompDtrace(Z, lambda = lambda, lambda.min.ratio = 0.01))
  CD_Omega_0_hat = CD_output$Theta[[which.min(CD_output$bic)]]
  CD_result = cbind(c(norm(CD_Omega_0_hat - Omega_0, "2"),
                      norm(CD_Omega_0_hat- Omega_0, "1"),
                      norm(CD_Omega_0_hat - Omega_0, "F"),
                      ROC_PR(CD_Omega_0_hat, Omega_0),CD_tt,min(eigen(CD_Omega_0_hat)$values)), CD_result)
  
  
  #----- gCoda -----#
  gCoda_tt<-ttime(gCoda_output <- gcoda(Z, nlambda = Nlambda, lambda.min.ratio = 0.01))
  gCoda_Omega_0_hat = gCoda_output$opt.icov
  gCoda_result = cbind(c(norm(gCoda_Omega_0_hat - Omega_0, "2"),
                         norm(gCoda_Omega_0_hat - Omega_0, "1"),
                         norm(gCoda_Omega_0_hat - Omega_0, "F"),
                         ROC_PR(gCoda_Omega_0_hat, Omega_0), gCoda_tt,min(eigen(gCoda_Omega_0_hat)$values)), gCoda_result)
  
  
  #----- Spiec_easi -----#
  Spiec_easi_tt<-ttime(Spiec_easi_output <- sparseiCov(Z, method = "glasso", lambda.min.ratio = 0.1))+ttime(Spiec_easi_Omega_hat <- as.matrix(huge.select(Spiec_easi_output, criterion = "stars", 
                                                                                                                                                        stars.thresh = 0.05, stars.subsample.ratio = 2 / 3, 
                                                                                                                                                        verbose = FALSE)$opt.icov))
  Spiec_easi_result = cbind(c(norm(Spiec_easi_Omega_hat - Omega_0, "2"),
                              norm(Spiec_easi_Omega_hat - Omega_0, "1"),
                              norm(Spiec_easi_Omega_hat - Omega_0, "F"),
                              ROC_PR(Spiec_easi_Omega_hat, Omega_0), Spiec_easi_tt,min(eigen(Spiec_easi_Omega_hat)$values)), Spiec_easi_result)
  
  cat(i, "Going on!\n")
}


#--------------------------------------- Summary -----------------------------------------------#
Care = round(cbind(apply(Care_result, 1, mean), apply(Care_result, 1, sd)), 5)
ada_bpm = round(cbind(apply(ada_bpm_result, 1, mean), apply(ada_bpm_result, 1, sd)), 5)
CD_trace = round(cbind(apply(CD_result, 1, mean), apply(CD_result, 1, sd)), 5)
gCoda = round(cbind(apply(gCoda_result, 1, mean), apply(gCoda_result, 1, sd)), 5)
Spiec_easi = round(cbind(apply(Spiec_easi_result, 1, mean), apply(Spiec_easi_result, 1, sd)), 5)


#--------------------------------------- Label names -------------------------------------------#
colnames(Care) = c("mean", "se"); colnames(ada_bpm) = c("mean", "se")
colnames(CD_trace) = c("mean", "se"); colnames(gCoda) = c("mean", "se")
colnames(Spiec_easi) = c("mean", "se")
rownames(Care) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(ada_bpm) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(CD_trace) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(gCoda) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(Spiec_easi) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")


#------------------------------------------ Save -----------------------------------------------#
Band_p_50 = list(Care = Care, ada_bpm = ada_bpm, CD_trace = CD_trace, 
                 gCoda = gCoda, Spiec_easi = Spiec_easi)
save(Band_p_50, file = "Band_p_50.RData")



#------------------------------ (Hub graph, model (b)) -------------------------------#
set.seed(12345)
n = 200
p = 50  # p = 50, 100, 200, 400，800
hub_am = make_graph("hub", p, p, enforce = TRUE, numHubs = p / 5)
nedge = sum(hub_am) / 2
temp = hub_am
edge_binary = rbinom(nedge, 1, 0.5)
edge_binary[edge_binary == 1] = 0.8; edge_binary[edge_binary == 0] = 0.5
temp[upper.tri(temp)] = 0; temp[temp != 0] = edge_binary
Omega_0 = temp + t(temp) + diag(runif(p, 1, 2))
Omega_0= (abs(min(eigen(Omega_0)$values)) + 0.01) * diag(p) + Omega_0
Sigma_0 = solve(Omega_0)
ratio = c(1/2, 7/12, 2/3, 3/4) 
Care_result = NULL; ada_bpm_result = NULL; CD_result = NULL
gCoda_result = NULL; Spiec_easi_result = NULL
for(i in 1 : 2){
  Y = rmvnorm(n, rep(0, p), Sigma_0)
  W = exp(Y)
  X = W / rowSums(W)
  G = diag(p) - matrix(1, p, p) / p
  Z = log(X) %*% G
  Nlambda = 50
  FsigmaF<-cov(Z)*(1-1/n);
  lam.max <- max(max(FsigmaF - diag(p)), -min(FsigmaF - diag(p)));
  lambda.min.ratio4=0.006
  lam.min4 <- lambda.min.ratio4 * lam.max;
  lam4 <- exp(seq(log(lam.max), log(lam.min4), length = Nlambda));
  
  
  
  #----- Care -----#
  Care_tt<-ttime(Care_output <- Care_est(X, nlambda = Nlambda, lambda_min = 0.01, ratio = 11/12, nsplit = 5))
  Care_Omega_hat = Care_output$Omega_hat
  Care_result = cbind(c(norm(Care_Omega_hat - Omega_0, "2"),
                        norm(Care_Omega_hat - Omega_0, "1"),
                        norm(Care_Omega_hat - Omega_0, "F"),
                        ROC_PR(Care_Omega_hat, Omega_0),Care_tt, min(eigen(Care_Omega_hat)$values)), Care_result)
  
  
  #----- ada_bpm -----#
  ada_bpm_tt<-ttime(ada_bpm_output <- ADA_BPM(Z,lambda =lam4, gamma=0.98))
  ada_bpm_Omega_hat=ada_bpm_output$Omega[which.min(ada_bpm_output$bic)][[1]]
  ada_bpm_Omega_hat=as.matrix(ada_bpm_Omega_hat)
  ada_bpm_result = cbind(c(Matrix::norm(ada_bpm_Omega_hat - Omega_0, "2"),
                        norm(ada_bpm_Omega_hat - Omega_0, "1"),
                        norm(ada_bpm_Omega_hat - Omega_0, "F"),
                        ROC_PR(ada_bpm_Omega_hat, Omega_0), ada_bpm_tt,min(eigen(ada_bpm_Omega_hat)$values)), ada_bpm_result)
  
  
  #----- CD-trace -----#
  lambda = seq(2 * (log(p) / n) ^ (1 / 2), 10 * (log(p) / n) ^ (1 / 2), length.out = Nlambda)
  CD_tt<-ttime(CD_output <- CompDtrace(Z, lambda = lambda, lambda.min.ratio = 0.01))
  CD_Omega_0_hat = CD_output$Theta[[which.min(CD_output$bic)]]
  CD_result = cbind(c(norm(CD_Omega_0_hat - Omega_0, "2"),
                      norm(CD_Omega_0_hat- Omega_0, "1"),
                      norm(CD_Omega_0_hat - Omega_0, "F"),
                      ROC_PR(CD_Omega_0_hat, Omega_0),CD_tt,min(eigen(CD_Omega_0_hat)$values)), CD_result)
  
  
  #----- gCoda -----#
  gCoda_tt<-ttime(gCoda_output <- gcoda(Z, nlambda = Nlambda, lambda.min.ratio = 0.01))
  gCoda_Omega_0_hat = gCoda_output$opt.icov
  gCoda_result = cbind(c(norm(gCoda_Omega_0_hat - Omega_0, "2"),
                         norm(gCoda_Omega_0_hat - Omega_0, "1"),
                         norm(gCoda_Omega_0_hat - Omega_0, "F"),
                         ROC_PR(gCoda_Omega_0_hat, Omega_0), gCoda_tt,min(eigen(gCoda_Omega_0_hat)$values)), gCoda_result)
  
  
  #----- Spiec_easi -----#
  Spiec_easi_tt<-ttime(Spiec_easi_output <- sparseiCov(Z, method = "glasso", lambda.min.ratio = 0.1))+ttime(Spiec_easi_Omega_hat <- as.matrix(huge.select(Spiec_easi_output, criterion = "stars", 
                                                                                                                                                          stars.thresh = 0.05, stars.subsample.ratio = 2 / 3, 
                                                                                                                                                          verbose = FALSE)$opt.icov))
  Spiec_easi_result = cbind(c(norm(Spiec_easi_Omega_hat - Omega_0, "2"),
                              norm(Spiec_easi_Omega_hat - Omega_0, "1"),
                              norm(Spiec_easi_Omega_hat - Omega_0, "F"),
                              ROC_PR(Spiec_easi_Omega_hat, Omega_0), Spiec_easi_tt,min(eigen(Spiec_easi_Omega_hat)$values)), Spiec_easi_result)
  
  cat(i, "Going on!\n")
}


#--------------------------------------- Summary -----------------------------------------------#
Care = round(cbind(apply(Care_result, 1, mean), apply(Care_result, 1, sd)), 5)
ada_bpm = round(cbind(apply(ada_bpm_result, 1, mean), apply(ada_bpm_result, 1, sd)), 5)
CD_trace = round(cbind(apply(CD_result, 1, mean), apply(CD_result, 1, sd)), 5)
gCoda = round(cbind(apply(gCoda_result, 1, mean), apply(gCoda_result, 1, sd)), 5)
Spiec_easi = round(cbind(apply(Spiec_easi_result, 1, mean), apply(Spiec_easi_result, 1, sd)), 5)


#--------------------------------------- Label names -------------------------------------------#
colnames(Care) = c("mean", "se"); colnames(ada_bpm) = c("mean", "se")
colnames(CD_trace) = c("mean", "se"); colnames(gCoda) = c("mean", "se")
colnames(Spiec_easi) = c("mean", "se")
rownames(Care) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(ada_bpm) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(CD_trace) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(gCoda) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(Spiec_easi) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")


#------------------------------------------ Save -----------------------------------------------#
Hub_p_50 = list(Care = Care, ada_bpm = ada_bpm, CD_trace = CD_trace, 
                gCoda = gCoda, Spiec_easi = Spiec_easi)
save(Hub_p_50, file = "Hub_p_50.RData")




#------------------------------ (Random graph, model (d)) -------------------------------#
set.seed(12345)
n = 200
p = 50  # p = 50, 100, 200, 400，800
prob = 4 / p 
random_graph = huge.generator(n, p, graph = "random", prob = prob, verbose = FALSE)
random_am = as.matrix(random_graph$theta)
nedge = sum(random_am) / 2; temp = random_am
edge_binary = rbinom(nedge, 1, 0.5)
edge_binary[edge_binary == 1] = 0.8; edge_binary[edge_binary == 0] = 0.5
temp[upper.tri(temp)] = 0; temp[temp != 0] = edge_binary
Omega_0 = temp + t(temp) + diag(runif(p, 1, 2))
Omega_0 = (abs(min(eigen(Omega_0)$values)) + 0.01) * diag(p) + Omega_0
Sigma_0 = solve(Omega_0)
ratio = c(1/2, 7/12, 2/3, 3/4) 
Care_result = NULL; ada_bpm_result = NULL; CD_result = NULL
gCoda_result = NULL; Spiec_easi_result = NULL
for(i in 1 : 2){
  Y = rmvnorm(n, rep(0, p), Sigma_0)
  W = exp(Y)
  X = W / rowSums(W)
  G = diag(p) - matrix(1, p, p) / p
  Z = log(X) %*% G
  Nlambda = 50
  FsigmaF<-cov(Z)*(1-1/n);
  lam.max <- max(max(FsigmaF - diag(p)), -min(FsigmaF - diag(p)));
  lambda.min.ratio4=0.001
  lam.min4 <- lambda.min.ratio4 * lam.max;
  lam4 <- exp(seq(log(lam.max), log(lam.min4), length = Nlambda));
  
  
  
  #----- Care -----#
  Care_tt<-ttime(Care_output <- Care_est(X, nlambda = Nlambda, lambda_min = 0.01, ratio = 11/12, nsplit = 5))
  Care_Omega_hat = Care_output$Omega_hat
  Care_result = cbind(c(norm(Care_Omega_hat - Omega_0, "2"),
                        norm(Care_Omega_hat - Omega_0, "1"),
                        norm(Care_Omega_hat - Omega_0, "F"),
                        ROC_PR(Care_Omega_hat, Omega_0),Care_tt, min(eigen(Care_Omega_hat)$values)), Care_result)
  
  
  #----- ada_bpm -----#
  ada_bpm_tt<-ttime(ada_bpm_output <- ADA_BPM(Z,lambda =lam4, gamma=0.999))
  ada_bpm_Omega_hat=ada_bpm_output$Omega[which.min(ada_bpm_output$bic)][[1]]
  ada_bpm_Omega_hat=as.matrix(ada_bpm_Omega_hat)
  ada_bpm_result = cbind(c(Matrix::norm(ada_bpm_Omega_hat - Omega_0, "2"),
                        norm(ada_bpm_Omega_hat - Omega_0, "1"),
                        norm(ada_bpm_Omega_hat - Omega_0, "F"),
                        ROC_PR(ada_bpm_Omega_hat, Omega_0), ada_bpm_tt,min(eigen(ada_bpm_Omega_hat)$values)), ada_bpm_result)
  
  
  #----- CD-trace -----#
  lambda = seq(2 * (log(p) / n) ^ (1 / 2), 10 * (log(p) / n) ^ (1 / 2), length.out = Nlambda)
  CD_tt<-ttime(CD_output <- CompDtrace(Z, lambda = lambda, lambda.min.ratio = 0.01))
  CD_Omega_0_hat = CD_output$Theta[[which.min(CD_output$bic)]]
  CD_result = cbind(c(norm(CD_Omega_0_hat - Omega_0, "2"),
                      norm(CD_Omega_0_hat- Omega_0, "1"),
                      norm(CD_Omega_0_hat - Omega_0, "F"),
                      ROC_PR(CD_Omega_0_hat, Omega_0),CD_tt,min(eigen(CD_Omega_0_hat)$values)), CD_result)
  
  
  #----- gCoda -----#
  gCoda_tt<-ttime(gCoda_output <- gcoda(Z, nlambda = Nlambda, lambda.min.ratio = 0.01))
  gCoda_Omega_0_hat = gCoda_output$opt.icov
  gCoda_result = cbind(c(norm(gCoda_Omega_0_hat - Omega_0, "2"),
                         norm(gCoda_Omega_0_hat - Omega_0, "1"),
                         norm(gCoda_Omega_0_hat - Omega_0, "F"),
                         ROC_PR(gCoda_Omega_0_hat, Omega_0), gCoda_tt,min(eigen(gCoda_Omega_0_hat)$values)), gCoda_result)
  
  
  #----- Spiec_easi -----#
  Spiec_easi_tt<-ttime(Spiec_easi_output <- sparseiCov(Z, method = "glasso", lambda.min.ratio = 0.1))+ttime(Spiec_easi_Omega_hat <- as.matrix(huge.select(Spiec_easi_output, criterion = "stars", 
                                                                                                                                                          stars.thresh = 0.05, stars.subsample.ratio = 2 / 3, 
                                                                                                                                                          verbose = FALSE)$opt.icov))
  Spiec_easi_result = cbind(c(norm(Spiec_easi_Omega_hat - Omega_0, "2"),
                              norm(Spiec_easi_Omega_hat - Omega_0, "1"),
                              norm(Spiec_easi_Omega_hat - Omega_0, "F"),
                              ROC_PR(Spiec_easi_Omega_hat, Omega_0), Spiec_easi_tt,min(eigen(Spiec_easi_Omega_hat)$values)), Spiec_easi_result)
  
  cat(i, "Going on!\n")
}


#--------------------------------------- Summary -----------------------------------------------#
Care = round(cbind(apply(Care_result, 1, mean), apply(Care_result, 1, sd)), 5)
ada_bpm = round(cbind(apply(ada_bpm_result, 1, mean), apply(ada_bpm_result, 1, sd)), 5)
CD_trace = round(cbind(apply(CD_result, 1, mean), apply(CD_result, 1, sd)), 5)
gCoda = round(cbind(apply(gCoda_result, 1, mean), apply(gCoda_result, 1, sd)), 5)
Spiec_easi = round(cbind(apply(Spiec_easi_result, 1, mean), apply(Spiec_easi_result, 1, sd)), 5)


#--------------------------------------- Label names -------------------------------------------#
colnames(Care) = c("mean", "se"); colnames(ada_bpm) = c("mean", "se")
colnames(CD_trace) = c("mean", "se"); colnames(gCoda) = c("mean", "se")
colnames(Spiec_easi) = c("mean", "se")
rownames(Care) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(ada_bpm) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(CD_trace) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(gCoda) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")
rownames(Spiec_easi) = c("S_norm", "L1_norm", "F_norm", "TPR", "FPR","time","min-eigen")


#------------------------------------------ Save -----------------------------------------------#
Random_p_50 = list(Care = Care, ada_bpm = ada_bpm, CD_trace = CD_trace, 
                  gCoda = gCoda, Spiec_easi = Spiec_easi)
save(Random_p_50, file = "Random_p_50.RData")




