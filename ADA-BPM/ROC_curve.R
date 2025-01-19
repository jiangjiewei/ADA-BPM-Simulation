#------------------------------------- Set path ---------------------------------------#
rm(list=ls())
setwd('~/ADA-BPM')


#----------------------------------- Load packages ------------------------------------#
library(ggplot2)
library(huge)
library(mvtnorm)
library(PRIMAL)
library(RColorBrewer)
library(MASS)
library('Rcpp')
library('Matrix')
#---------------------------------- Import functions ----------------------------------#
source("Care.R")
source("Oracle.R")
source("CD-trace.R")
source("gCoda.R")
source("SPIEC-EASI.R")
source("LNM.R")
source("Naive_methods.R")
source("Summary_roc.R")
sourceCpp("Trf.cpp")
sourceCpp('adabpm.cpp')
source("ADA-BPM.R")

#------------------------- ROC curves for compositional data --------------------------#
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#

#------------------------------ (Band graph, model (a)) -------------------------------#
set.seed(12345)
n = 200
p = 50  # p = 50, 100, 200, 400
band_graph = huge.generator(n, p, graph = "band", g = 2, verbose = FALSE)
band_am = as.matrix(band_graph$theta)
band_graph1 = huge.generator(n, p, graph = "band", g = 1, verbose = FALSE)
band_am1 = as.matrix(band_graph1$theta)
Omega_0 = band_am * 0.5 + band_am1 * 0.3 + diag(runif(p, 1, 2))
Omega_0 = (abs(min(eigen(Omega_0)$values)) + 0.01) * diag(p) + Omega_0
Sigma_0 = solve(Omega_0)



Care_roc = list(); face_roc = list(); CD_roc = list();
gCoda_roc = list(); Spieceasi_roc = list();

for(i in 1 : 10){
  set.seed(i)
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
  Care_path = Care_sp(X, nlambda = 50, lambda_min = 0.001)
  Care_roc[[i]] = summary_ROC_PR(Care_path, Omega_0)
  
  
  #----- face -----#
  face_output <- ADA_BPM(Z,lambda =lam4, gamma=0.5)
  face_roc[[i]] = summary_ROC_PR(face_output$path, Omega_0)
  
  
  
  #----- CD-trace -----#
  CD_output =  CompDtrace(Z, nlambda = 50, rho = 5, lambda.min.ratio = 0.001)
  CD_roc[[i]] = summary_ROC_PR(CD_output$path, Omega_0)
  
  
  #----- gCoda -----#
  gCoda_output = gcoda(Z, nlambda = 50, lambda.min.ratio = 0.001)
  gCoda_roc[[i]] = summary_ROC_PR(gCoda_output$path, Omega_0)
  
  
  #----- Spiec_easi -----#
  Spieceasi_output = sparseiCov(Z, method = "glasso", 
                                nlambda = 50, lambda.min.ratio = 0.001)
  Spieceasi_roc[[i]] = summary_ROC_PR(Spieceasi_output$path, Omega_0)
  
  
  cat(i, "Going on!\n")
}

#--------------------------------------- Summary -----------------------------------------------#
Care_TPR = c(rowMeans(sapply(Care_roc, function(y)y$TPR)), 1.0)
Care_FPR = c(rowMeans(sapply(Care_roc, function(y)y$FPR)), 1.0)
face_TPR = c(rowMeans(sapply(face_roc, function(y)y$TPR)), 1.0)
face_FPR = c(rowMeans(sapply(face_roc, function(y)y$FPR)), 1.0)
CD_TPR = c(rowMeans(sapply(CD_roc, function(y)y$TPR)), 1.0)
CD_FPR = c(rowMeans(sapply(CD_roc, function(y)y$FPR)), 1.0)
gCoda_TPR = c(rowMeans(sapply(gCoda_roc, function(y)y$TPR)), 1.0)
gCoda_FPR = c(rowMeans(sapply(gCoda_roc, function(y)y$FPR)), 1.0)
Spieceasi_TPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$TPR)), 1.0)
Spieceasi_FPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$FPR)), 1.0)
Method = rep(c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI"), c(50, 50, 50, 50, 50) + 1)
Method = factor(Method)
TPR = c(Care_TPR, face_TPR, CD_TPR, gCoda_TPR, Spieceasi_TPR)
FPR = c(Care_FPR, face_FPR, CD_FPR, gCoda_FPR, Spieceasi_FPR)
Band_ROC = data.frame(cbind(TPR, FPR))
Band_ROC$Method = Method

#---------------------------------------- Plot for ROC ----------------------------------------#
theme_set(theme_bw())
Band_figure = ggplot(Band_ROC, aes(x = FPR, y = TPR)) + 
  geom_line(linewidth = 0.7, aes(color = Method, linetype = Method)) + 
  ggtitle(expression(paste("Band graph (", italic("p "), "= 50)"))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "twodash", "dotted"), 
                        breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_color_manual(values = c("#E31A1C", "#FF7F00", "#1F78B4", "#33A02C", "#FB9A99"),
                     breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_x_continuous(name = "False positive rate", breaks = seq(0.0, 1.0, 0.2)) + 
  scale_y_continuous(name = "True positive rate", breaks = seq(0.0, 1.0, 0.2))
Band_figure + 
  theme(axis.text.x= element_text(color = "black"), 
        axis.text.y= element_text(color = "black")) +
  theme(text = element_text(size = 20), plot.margin = margin(0, 0.3, 0.2, 0.2, "cm"),
        plot.title = element_text(size = 20, vjust = 0.1, hjust = 0.5)) + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0),
        legend.box.margin=margin(c(10, 10, 10, 10))) +
  theme(panel.grid=element_blank())


#------------------------------ (Hub graph, model (b)) -------------------------------#
set.seed(12345)
n = 200
p = 50  # p = 50, 100, 200, 400
hub_am = make_graph("hub", p, p, enforce = TRUE, numHubs = p / 5)
nedge = sum(hub_am) / 2
temp = hub_am
edge_binary = rbinom(nedge, 1, 0.5)
edge_binary[edge_binary == 1] = 0.8; edge_binary[edge_binary == 0] = 0.5
temp[upper.tri(temp)] = 0; temp[temp != 0] = edge_binary
Omega_0 = temp + t(temp) + diag(runif(p, 1, 2))
Omega_0= (abs(min(eigen(Omega_0)$values)) + 0.01) * diag(p) + Omega_0
Sigma_0 = solve(Omega_0)
Care_roc = list(); face_roc = list(); CD_roc = list()
gCoda_roc = list(); Spieceasi_roc = list()
for(i in 1 : 10){
  set.seed(i)
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
  Care_path = Care_sp(X, nlambda = 50, lambda_min = 0.001)
  Care_roc[[i]] = summary_ROC_PR(Care_path, Omega_0)
  
  
  #----- face -----#
  face_output <- ADA_BPM(Z,lambda =lam4, gamma=0.5)
  face_roc[[i]] = summary_ROC_PR(face_output$path, Omega_0)
  
  
  #----- CD-trace -----#
  CD_output =  CompDtrace(Z, nlambda = 50, rho = 5, lambda.min.ratio = 0.001)
  CD_roc[[i]] = summary_ROC_PR(CD_output$path, Omega_0)
  
  
  #----- gCoda -----#
  gCoda_output = gcoda(Z, nlambda = 50, lambda.min.ratio = 0.001)
  gCoda_roc[[i]] = summary_ROC_PR(gCoda_output$path, Omega_0)
  
  
  #----- Spiec_easi -----#
  Spieceasi_output = sparseiCov(Z, method = "glasso", 
                                nlambda = 50, lambda.min.ratio = 0.001)
  Spieceasi_roc[[i]] = summary_ROC_PR(Spieceasi_output$path, Omega_0)
  
  
  cat(i, "Going on!\n")
}

#--------------------------------------- Summary -----------------------------------------------#
Care_TPR = c(rowMeans(sapply(Care_roc, function(y)y$TPR)), 1.0)
Care_FPR = c(rowMeans(sapply(Care_roc, function(y)y$FPR)), 1.0)
face_TPR = c(rowMeans(sapply(face_roc, function(y)y$TPR)), 1.0)
face_FPR = c(rowMeans(sapply(face_roc, function(y)y$FPR)), 1.0)
CD_TPR = c(rowMeans(sapply(CD_roc, function(y)y$TPR)), 1.0)
CD_FPR = c(rowMeans(sapply(CD_roc, function(y)y$FPR)), 1.0)
gCoda_TPR = c(rowMeans(sapply(gCoda_roc, function(y)y$TPR)), 1.0)
gCoda_FPR = c(rowMeans(sapply(gCoda_roc, function(y)y$FPR)), 1.0)
Spieceasi_TPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$TPR)), 1.0)
Spieceasi_FPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$FPR)), 1.0)
Method = rep(c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI"), c(50, 50, 50, 50, 50) + 1)
Method = factor(Method)
TPR = c(Care_TPR, face_TPR, CD_TPR, gCoda_TPR, Spieceasi_TPR)
FPR = c(Care_FPR, face_FPR, CD_FPR, gCoda_FPR, Spieceasi_FPR)
Hub_ROC = data.frame(cbind(TPR, FPR))
Hub_ROC$Method = Method

#---------------------------------------- Plot for ROC ----------------------------------------#
theme_set(theme_bw())
Hub_figure = ggplot(Hub_ROC, aes(x = FPR, y = TPR)) + 
  geom_line(linewidth = 0.7, aes(color = Method, linetype = Method)) + 
  ggtitle(expression(paste("Hub graph (", italic("p "), "= 50)"))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "twodash", "dotted"), 
                        breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_color_manual(values = c("#E31A1C", "#FF7F00", "#1F78B4", "#33A02C", "#FB9A99"),
                     breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_x_continuous(name = "False positive rate", breaks = seq(0.0, 1.0, 0.2)) + 
  scale_y_continuous(name = "True positive rate", breaks = seq(0.0, 1.0, 0.2))
Hub_figure + 
  theme(axis.text.x= element_text(color = "black"), 
        axis.text.y= element_text(color = "black")) +
  theme(text = element_text(size = 20), plot.margin = margin(0, 0.3, 0.2, 0.2, "cm"),
        plot.title = element_text(size = 20, vjust = 0.1, hjust = 0.5)) + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0),
        legend.box.margin=margin(c(10, 10, 10, 10))) +
  theme(panel.grid=element_blank())


#------------------------------ (Block graph, model (c)) -------------------------------#
set.seed(12345)
n = 200
p = 50  # p = 50, 100, 200, 400
nblock = 5 
prob = 20 / p 
block_graph = huge.generator(n, p, graph = "cluster", 
                             g = nblock, prob = prob, verbose = FALSE)
block_am = as.matrix(block_graph$theta)
nedge = sum(block_am) / 2; temp = block_am
edge_binary = rbinom(nedge, 1, 0.5)
edge_binary[edge_binary == 1] = 0.8; edge_binary[edge_binary == 0] = 0.5
temp[upper.tri(temp)] = 0; temp[temp != 0] = edge_binary
Omega_0 = temp + t(temp) + diag(runif(p, 1, 2))
Omega_0 = (abs(min(eigen(Omega_0)$values)) + 0.01) * diag(p) + Omega_0
Sigma_0 = solve(Omega_0)
Care_roc = list(); face_roc = list(); CD_roc = list()
gCoda_roc = list(); Spieceasi_roc = list()
for(i in 1 : 10){
  set.seed(i)
  Y = rmvnorm(n, rep(0, p), Sigma_0)
  W = exp(Y)
  X = W / rowSums(W)
  G = diag(p) - matrix(1, p, p) / p
  Z = log(X) %*% G
  
  
  #----- Care -----#
  Care_path = Care_sp(X, nlambda = 50, lambda_min = 0.001)
  Care_roc[[i]] = summary_ROC_PR(Care_path, Omega_0)
  
  
  #----- face -----#
  face_output <- ADA_BPM(Z, gamma=0.5)
  face_roc[[i]] = summary_ROC_PR(face_output$path, Omega_0)
  
  
  #----- CD-trace -----#
  CD_output =  CompDtrace(Z, nlambda = 50, rho = 5, lambda.min.ratio = 0.001)
  CD_roc[[i]] = summary_ROC_PR(CD_output$path, Omega_0)
  
  
  #----- gCoda -----#
  gCoda_output = gcoda(Z, nlambda = 50, lambda.min.ratio = 0.001)
  gCoda_roc[[i]] = summary_ROC_PR(gCoda_output$path, Omega_0)
  
  
  #----- Spiec_easi -----#
  Spieceasi_output = sparseiCov(Z, method = "glasso", 
                                nlambda = 50, lambda.min.ratio = 0.001)
  Spieceasi_roc[[i]] = summary_ROC_PR(Spieceasi_output$path, Omega_0)
  
  
  cat(i, "Going on!\n")
}

#--------------------------------------- Summary -----------------------------------------------#
Care_TPR = c(rowMeans(sapply(Care_roc, function(y)y$TPR)), 1.0)
Care_FPR = c(rowMeans(sapply(Care_roc, function(y)y$FPR)), 1.0)
face_TPR = c(rowMeans(sapply(face_roc, function(y)y$TPR)), 1.0)
face_FPR = c(rowMeans(sapply(face_roc, function(y)y$FPR)), 1.0)
CD_TPR = c(rowMeans(sapply(CD_roc, function(y)y$TPR)), 1.0)
CD_FPR = c(rowMeans(sapply(CD_roc, function(y)y$FPR)), 1.0)
gCoda_TPR = c(rowMeans(sapply(gCoda_roc, function(y)y$TPR)), 1.0)
gCoda_FPR = c(rowMeans(sapply(gCoda_roc, function(y)y$FPR)), 1.0)
Spieceasi_TPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$TPR)), 1.0)
Spieceasi_FPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$FPR)), 1.0)
Method = rep(c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI"), c(50, 50, 50, 50, 50) + 1)
Method = factor(Method)
TPR = c(Care_TPR, face_TPR, CD_TPR, gCoda_TPR, Spieceasi_TPR)
FPR = c(Care_FPR, face_FPR, CD_FPR, gCoda_FPR, Spieceasi_FPR)
Block_ROC = data.frame(cbind(TPR, FPR))
Block_ROC$Method = Method

#---------------------------------------- Plot for ROC ----------------------------------------#
theme_set(theme_bw())
Block_figure = ggplot(Block_ROC, aes(x = FPR, y = TPR)) + 
  geom_line(linewidth = 0.7, aes(color = Method, linetype = Method)) + 
  ggtitle(expression(paste("Block graph (", italic("p "), "= 50)"))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "twodash", "dotted"), 
                        breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_color_manual(values = c("#E31A1C", "#FF7F00", "#1F78B4", "#33A02C", "#FB9A99"),
                     breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_x_continuous(name = "False positive rate", breaks = seq(0.0, 1.0, 0.2)) + 
  scale_y_continuous(name = "True positive rate", breaks = seq(0.0, 1.0, 0.2))
Block_figure + 
  theme(axis.text.x= element_text(color = "black"), 
        axis.text.y= element_text(color = "black")) +
  theme(text = element_text(size = 20), plot.margin = margin(0, 0.3, 0.2, 0.2, "cm"),
        plot.title = element_text(size = 20, vjust = 0.1, hjust = 0.5)) + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0),
        legend.box.margin=margin(c(10, 10, 10, 10))) +
  theme(panel.grid=element_blank())


#------------------------------ (Random graph, model (d)) -------------------------------#
set.seed(12345)
n = 200
p = 50  # p = 50, 100, 200, 400
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
Care_roc = list(); face_roc = list(); CD_roc = list()
gCoda_roc = list(); Spieceasi_roc = list()
for(i in 1 : 10){
  set.seed(i)
  Y = rmvnorm(n, rep(0, p), Sigma_0)
  W = exp(Y)
  X = W / rowSums(W)
  G = diag(p) - matrix(1, p, p) / p
  Z = log(X) %*% G
  
  
  #----- Care -----#
  Care_path = Care_sp(X, nlambda = 50, lambda_min = 0.001)
  Care_roc[[i]] = summary_ROC_PR(Care_path, Omega_0)
  
  
  #----- face -----#
  face_output <- ADA_BPM(Z,lambda =lam4, gamma=0.5)
  face_roc[[i]] = summary_ROC_PR(face_output$path, Omega_0)
  
  
  #----- CD-trace -----#
  CD_output =  CompDtrace(Z, nlambda = 50, rho = 5, lambda.min.ratio = 0.001)
  CD_roc[[i]] = summary_ROC_PR(CD_output$path, Omega_0)
  
  
  #----- gCoda -----#
  gCoda_output = gcoda(Z, nlambda = 50, lambda.min.ratio = 0.001)
  gCoda_roc[[i]] = summary_ROC_PR(gCoda_output$path, Omega_0)
  
  
  #----- Spiec_easi -----#
  Spieceasi_output = sparseiCov(Z, method = "glasso", 
                                nlambda = 50, lambda.min.ratio = 0.001)
  Spieceasi_roc[[i]] = summary_ROC_PR(Spieceasi_output$path, Omega_0)
  
  
  cat(i, "Going on!\n")
}

#--------------------------------------- Summary -----------------------------------------------#
Care_TPR = c(rowMeans(sapply(Care_roc, function(y)y$TPR)), 1.0)
Care_FPR = c(rowMeans(sapply(Care_roc, function(y)y$FPR)), 1.0)
face_TPR = c(rowMeans(sapply(face_roc, function(y)y$TPR)), 1.0)
face_FPR = c(rowMeans(sapply(face_roc, function(y)y$FPR)), 1.0)
CD_TPR = c(rowMeans(sapply(CD_roc, function(y)y$TPR)), 1.0)
CD_FPR = c(rowMeans(sapply(CD_roc, function(y)y$FPR)), 1.0)
gCoda_TPR = c(rowMeans(sapply(gCoda_roc, function(y)y$TPR)), 1.0)
gCoda_FPR = c(rowMeans(sapply(gCoda_roc, function(y)y$FPR)), 1.0)
Spieceasi_TPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$TPR)), 1.0)
Spieceasi_FPR = c(rowMeans(sapply(Spieceasi_roc, function(y)y$FPR)), 1.0)
Method = rep(c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI"), c(50, 50, 50, 50, 50) + 1)
Method = factor(Method)
TPR = c(Care_TPR, face_TPR, CD_TPR, gCoda_TPR, Spieceasi_TPR)
FPR = c(Care_FPR, face_FPR, CD_FPR, gCoda_FPR, Spieceasi_FPR)
Random_ROC = data.frame(cbind(TPR, FPR))
Random_ROC$Method = Method

#---------------------------------------- Plot for ROC ----------------------------------------#
theme_set(theme_bw())
Random_figure = ggplot(Random_ROC, aes(x = FPR, y = TPR)) + 
  geom_line(linewidth = 0.7, aes(color = Method, linetype = Method)) + 
  ggtitle(expression(paste("Random graph (", italic("p "), "= 50)"))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "twodash", "dotted"), 
                        breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_color_manual(values = c("#E31A1C", "#FF7F00", "#1F78B4", "#33A02C", "#FB9A99"),
                     breaks = c("CARE", "face", "CD-trace", "gCoda", "SPIEC-EASI")) + 
  scale_x_continuous(name = "False positive rate", breaks = seq(0.0, 1.0, 0.2)) + 
  scale_y_continuous(name = "True positive rate", breaks = seq(0.0, 1.0, 0.2))
Random_figure + 
  theme(axis.text.x= element_text(color = "black"), 
        axis.text.y= element_text(color = "black")) +
  theme(text = element_text(size = 20), plot.margin = margin(0, 0.3, 0.2, 0.2, "cm"),
        plot.title = element_text(size = 20, vjust = 0.1, hjust = 0.5)) + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0),
        legend.box.margin=margin(c(10, 10, 10, 10))) +
  theme(panel.grid=element_blank())

