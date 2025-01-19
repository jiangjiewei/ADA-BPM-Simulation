#------------------------------------- Set path ---------------------------------------#
rm(list=ls())
setwd('~/ada-bpm-real-data')

#----------------------------------- Load packages ------------------------------------#
library(ggplot2)
library(huge)
library(mvtnorm)
library(PRIMAL)
library(RColorBrewer)
library(scales)
library(igraph)
library(network)
library('MASS')
library('Rcpp')
library('Matrix')

#---------------------------------- Import functions ----------------------------------#
source("Care.R")
source("gCoda.R")
source("CD-trace.R")
source("SPIEC-EASI.R")
source("Summary_network.R")
source("basic_functions.R")
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(corpcor)
source("tune_proximal_gradient.R")
source("proximal_gradient.R")
sourceCpp("Trf.cpp")
sourceCpp('adabpm.cpp')
source("ada-bpm.R")

#------------------------------------ Import data -------------------------------------#
load("Gut Microbiome.RData")


#------------------------- Partition of lean and obese groups -------------------------#
names = row.names(countdata)
lean_index = NULL; obese_index = NULL
for(i in 1 : length(names)){
  if(names[i] %in% row.names(lean)){
    lean_index = c(lean_index, i)
  }
  if(names[i] %in% row.names(obese)){
    obese_index = c(obese_index, i)
  }
}


#--------------------------- Screen target genus and phylum --------------------------#
lean_num_appear = rep(0, ncol(countdata)); obese_num_appear = rep(0, ncol(countdata))
for(i in 1 : ncol(countdata)){
  lean_num_appear[i] = length(which(lean[, i] != 0))
  obese_num_appear[i] = length(which(obese[, i] != 0))
}
index = ifelse(lean_num_appear >= 4, 1, 0) * ifelse(obese_num_appear >= 4, 1, 0)
(p = length(which(index == 1)))
new_countdata = countdata[, which(index == 1)]


#------------------------------ Construct nodes dataset ------------------------------#
splitname = sapply(colnames(new_countdata), strsplit, "[.]")
genus = rep(0, p); phylum = rep(0, p)
for (i in 1 : p){
  genus[i] = splitname[[i]][length(splitname[[i]])]
  phylum[i] = splitname[[i]][2]
}
phylum_type = rep(0, p)
phylum_type[phylum == "Actinobacteria"] = 1
phylum_type[phylum == "Bacteroidetes"] = 2
phylum_type[phylum == "Firmicutes"] = 3
phylum_type[phylum == "Proteobacteria"] = 4
Nodes = data.frame(cbind(genus, phylum, phylum_type))


#-------------------- Plot the singular values of MLE (X_hat) -----------------------#
composition = new_countdata / rowSums(new_countdata)
values = svd(composition)$d
id = 1 : p
singular_values = data.frame(cbind(id, values))
theme_set(theme_bw())
ggplot(singular_values, aes(x = id, y = values)) + geom_line(linewidth = 2) +
  scale_x_continuous(name = "Number of singular values") +
  scale_y_continuous(name = "Singular values", breaks = seq(0, 6, 1)) +
  theme(axis.text.x= element_text(color = "black"), 
        axis.text.y= element_text(color = "black")) + 
  theme(text = element_text(size = 15)) + theme(panel.grid=element_blank())


#------------------- Composition estimation by Cao et al. (2020) --------------------#
# PG_res = autoTuneProxGradient(W = new_countdata, n_grid = 5)
# composition_hat = PG_res$X_hat
# save(composition_hat, file = "composition_hat.RData")

load("composition_hat.RData")


#------------------------ CLR transformation of compositions ------------------------#
lean_composition = composition_hat[lean_index, ]
obese_composition = composition_hat[obese_index, ]
G = diag(p) - matrix(1, p, p) / p
lean_composition_clr = log(lean_composition) %*% G
obese_composition_clr = log(obese_composition) %*% G
Nlambda = 50

n_lean<-nrow(lean_composition_clr)
FsigmaF1<-cov(lean_composition_clr)*(1-1/n_lean);
lam.max <- max(max(FsigmaF1 - diag(p)), -min(FsigmaF1 - diag(p)));
lambda.min.ratio4=0.01
lam.min4 <- lambda.min.ratio4 * lam.max;
lam4 <- exp(seq(log(lam.max), log(lam.min4), length = Nlambda));

n_obese<-nrow(obese_composition_clr)
FsigmaF2<-cov(lean_composition_clr)*(1-1/n_obese);
lam.max1 <- max(max(FsigmaF2 - diag(p)), -min(FsigmaF2 - diag(p)));
lambda.min.ratio5=0.011
lam.min5 <- lambda.min.ratio5 * lam.max1;
lam5 <- exp(seq(log(lam.max1), log(lam.min5), length =Nlambda ));

#-------------------------------------- CARE ----------------------------------------#
set.seed(123456)

Care_lean_output = Care_est(lean_composition, nlambda = Nlambda, 
                            lambda_min = 0.01, ratio = 0.75, nsplit = 10)
Care_lean_Omega_0_hat = Care_lean_output$Omega_hat

Care_obese_output = Care_est(obese_composition, nlambda = Nlambda, 
                             lambda_min = 0.01, ratio = 0.75, nsplit = 10)
Care_obese_Omega_0_hat = Care_obese_output$Omega_hat

# compute the numbers of positive and negative edges
(Care_lean_Num_strength = Strength(Care_lean_Omega_0_hat))
# 9 + and 19 -
(Care_obese_Num_strength = Strength(Care_obese_Omega_0_hat))
# 4 + and 8 -

#------------------------------------ CD-trace --------------------------------------#
set.seed(123456)

CD_lean_output = CompDtrace(lean_composition_clr, nlambda = Nlambda, 
                            lambda.min.ratio = 0.01)
CD_obese_output = CompDtrace(obese_composition_clr, nlambda = Nlambda, 
                             lambda.min.ratio = 0.01)

CD_lean_Omega_0_hat = (CD_lean_output$Theta)[[which.min(CD_lean_output$bic)]]
CD_obese_Omega_0_hat = (CD_obese_output$Theta)[[which.min(CD_obese_output$bic)]]

# compute the numbers of positive and negative edges
(CD_lean_Num_strength = Strength(CD_lean_Omega_0_hat))
# 9 + and 10 -
(CD_obese_Num_strength = Strength(CD_obese_Omega_0_hat))
# 16 + and 15 -



#------------------------------------ ada-bpm --------------------------------------#
set.seed(123456)

face_lean_output <- ADA_BPM(lean_composition_clr,lambda =lam4, gamma=0.999)
face_obese_output <- ADA_BPM(obese_composition_clr,lambda =lam5, gamma=0.999)
face_lean_Omega_hat=face_lean_output$Omega[[which.min(face_lean_output$bic)]]
face_obese_Omega_hat=face_obese_output$Omega[[which.min(face_obese_output$bic)]]

# compute the numbers of positive and negative edges
(face_lean_Num_strength = Strength(face_lean_Omega_hat))
# 6 + and 13 -
(face_obese_Num_strength = Strength(face_obese_Omega_hat))
# 7 + and 12 -





#------------------------------------ gCoda ----------------------------------------#
set.seed(123456)

gCoda_lean_Omega_0_hat = gcoda(lean_composition_clr, lambda.min.ratio = 0.01, 
                               nlambda = Nlambda)$opt.icov
gCoda_obese_Omega_0_hat = gcoda(obese_composition_clr, lambda.min.ratio = 0.01, 
                                nlambda = Nlambda)$opt.icov

# compute the numbers of positive and negative edges
(gCoda_lean_Num_strength = Strength(gCoda_lean_Omega_0_hat))
# 12 + and 16 -
(gCoda_obese_Num_strength = Strength(gCoda_obese_Omega_0_hat))
# 4 + and 6 -


#--------------------------------- SPIEC-EASI -------------------------------------#
set.seed(123456)

SE_gl_lean_output = sparseiCov(lean_composition_clr, method = "glasso", 
                               lambda.min.ratio = 0.01, nlambda = Nlambda)
SE_gl_obese_output = sparseiCov(obese_composition_clr, method = "glasso", 
                                lambda.min.ratio = 0.01, nlambda = Nlambda)
SE_gl_lean_Omega_0_hat = as.matrix(huge.select(SE_gl_lean_output, 
                                               criterion = "stars", 
                                               verbose = FALSE)$opt.icov)
SE_gl_obese_Omega_0_hat = as.matrix(huge.select(SE_gl_obese_output, 
                                                criterion = "stars", 
                                                verbose = FALSE)$opt.icov)

# compute the numbers of positive and negative edges
(SE_gl_lean_Num_strength = Strength(SE_gl_lean_Omega_0_hat))
# 12 + and 18 -
(SE_gl_obese_Num_strength = Strength(SE_gl_obese_Omega_0_hat))
# 13 + and 9 -


#------------------ Compute network stability and edge replicates -----------------#
lean_list = list(Care_lean_Omega_0_hat, CD_lean_Omega_0_hat, face_lean_Omega_hat,
                 gCoda_lean_Omega_0_hat, SE_gl_lean_Omega_0_hat)
obese_list = list(Care_obese_Omega_0_hat, CD_obese_Omega_0_hat, face_obese_Omega_hat,
                  gCoda_obese_Omega_0_hat, SE_gl_obese_Omega_0_hat)


#lean_summary = Subsampling(lean_list, lean_composition, pro = 0.8, num = 50)
obese_summary = Subsampling(obese_list, obese_composition, pro = 0.8, num = 20)

#save(lean_summary, file = "lean_summary.RData")
#save(obese_summary, file = "obese_summary.RData")

load("lean_summary.RData")
load("obese_summary.RData")

#load("networks.RData")





#------------------------------------- CARE ---------------------------------------#
Care_lean_Edges = Edges_data(Care_lean_Omega_0_hat, Nodes)
Care_lean_Edges$rep = lean_summary[[2]]
#Care_lean_Edges$rep = networks$lean$Care_edge_rep
Care_lean_Edges = Care_lean_Edges[which(Care_lean_Edges$rep >= 40), ]
row.names(Care_lean_Edges) = NULL
(Care_lean_positive = length(which(Care_lean_Edges$sign == 1)))  # 2+
(Care_lean_negative = length(which(Care_lean_Edges$sign == -1))) # 5 -
# View(Care_lean_Edges)

Care_obese_Edges = Edges_data(Care_obese_Omega_0_hat, Nodes)
Care_obese_Edges$rep = obese_summary[[2]]
#Care_obese_Edges$rep = networks$obese$Care_edge_rep
Care_obese_Edges = Care_obese_Edges[which(Care_obese_Edges$rep >= 16), ]
row.names(Care_obese_Edges) = NULL
(Care_obese_positive = length(which(Care_obese_Edges$sign == 1)))  # 0 +
(Care_obese_negative = length(which(Care_obese_Edges$sign == -1))) # 1 -
# View(Care_obese_Edges)


#---------------------------------- CD-trace ---------------------------------------#
CD_lean_Edges = Edges_data(CD_lean_Omega_0_hat, Nodes)
CD_lean_Edges$rep = lean_summary[[3]]
#CD_lean_Edges$rep = networks$lean$CD_trace_edge_rep
CD_lean_Edges = CD_lean_Edges[which(CD_lean_Edges$rep >= 40), ]
row.names(CD_lean_Edges) = NULL
(CD_lean_positive = length(which(CD_lean_Edges$sign == 1)))  # 0 +
(CD_lean_negative = length(which(CD_lean_Edges$sign == -1))) # 1 -
# View(CD_lean_Edges)

CD_obese_Edges = Edges_data(CD_obese_Omega_0_hat, Nodes)
CD_obese_Edges$rep = obese_summary[[3]]
#CD_obese_Edges$rep = networks$obese$CD_trace_edge_rep
CD_obese_Edges = CD_obese_Edges[which(CD_obese_Edges$rep >= 16), ]
row.names(CD_obese_Edges) = NULL
(CD_obese_positive = length(which(CD_obese_Edges$sign == 1)))  # 0 +
(CD_obese_negative = length(which(CD_obese_Edges$sign == -1))) # 1 -
# View(CD_obese_Edges)


#---------------------------------- ada-bpm ---------------------------------------#
face_lean_Edges = Edges_data(face_lean_Omega_hat, Nodes)
face_lean_Edges$rep = lean_summary[[4]]
face_lean_Edges = face_lean_Edges[which(face_lean_Edges$rep >= 40), ]
row.names(face_lean_Edges) = NULL
(face_lean_positive = length(which(face_lean_Edges$sign == 1)))  #5 +
(face_lean_negative = length(which(face_lean_Edges$sign == -1))) #10 -
#View(face_lean_Edges)

face_obese_Edges = Edges_data(face_obese_Omega_hat, Nodes)
face_obese_Edges$rep = obese_summary[[4]]
face_obese_Edges = face_obese_Edges[which(face_obese_Edges$rep >= 16), ]
row.names(face_obese_Edges) = NULL
(face_obese_positive = length(which(face_obese_Edges$sign == 1)))  # 3 +
(face_obese_negative = length(which(face_obese_Edges$sign == -1))) # 8 -
#View(face_obese_Edges)

#------------------------------------ gCoda ----------------------------------------#
gCoda_lean_Edges = Edges_data(gCoda_lean_Omega_0_hat, Nodes)
gCoda_lean_Edges$rep = lean_summary[[5]]
#gCoda_lean_Edges$rep = networks$lean$gCoda_edge_rep
gCoda_lean_Edges = gCoda_lean_Edges[which(gCoda_lean_Edges$rep >= 40), ]
row.names(gCoda_lean_Edges) = NULL
(gCoda_lean_positive = length(which(gCoda_lean_Edges$sign == 1)))  # 4 +
(gCoda_lean_negative = length(which(gCoda_lean_Edges$sign == -1))) # 7 -
# View(gCoda_lean_Edges)

gCoda_obese_Edges = Edges_data(gCoda_obese_Omega_0_hat, Nodes)
gCoda_obese_Edges$rep = obese_summary[[5]]
#gCoda_obese_Edges$rep = networks$obese$gCoda_edge_rep
gCoda_obese_Edges = gCoda_obese_Edges[which(gCoda_obese_Edges$rep >= 16), ]
row.names(gCoda_obese_Edges) = NULL
(gCoda_obese_positive = length(which(gCoda_obese_Edges$sign == 1)))  # 0 +
(gCoda_obese_negative = length(which(gCoda_obese_Edges$sign == -1))) # 2 -
# View(gCoda_obese_Edges)


#---------------------------------- SPIEC-EASI -------------------------------------#
SE_gl_lean_Edges = Edges_data(SE_gl_lean_Omega_0_hat, Nodes)
SE_gl_lean_Edges$rep = lean_summary[[6]]
#SE_gl_lean_Edges$rep = networks$lean$Spiec_easi_edge_rep
SE_gl_lean_Edges = SE_gl_lean_Edges[which(SE_gl_lean_Edges$rep >= 40), ]
row.names(SE_gl_lean_Edges) = NULL
(SE_gl_lean_positive = length(which(SE_gl_lean_Edges$sign == 1)))  # 4 +
(SE_gl_lean_negative = length(which(SE_gl_lean_Edges$sign == -1))) # 5 -
# View(SE_gl_lean_Edges)

SE_gl_obese_Edges = Edges_data(SE_gl_obese_Omega_0_hat, Nodes)
SE_gl_obese_Edges$rep = obese_summary[[6]]
#SE_gl_obese_Edges$rep = networks$obese$Spiec_easi_edge_rep
SE_gl_obese_Edges = SE_gl_obese_Edges[which(SE_gl_obese_Edges$rep >= 20), ]
row.names(SE_gl_obese_Edges) = NULL
(SE_gl_obese_positive = length(which(SE_gl_obese_Edges$sign == 1)))  # 1 +
(SE_gl_obese_negative = length(which(SE_gl_obese_Edges$sign == -1))) # 2 -
# View(SE_gl_obese_Edges)






#----------------- Sort genera by the relative abundance and phylum -----------------#
re_abundance = colMeans(composition_hat)
re_abundance = as.numeric(round(log(re_abundance * 100000)))
Nodes$re_abundance = re_abundance

Nodes1 = Nodes[which(Nodes$phylum_type == "1"), ]
Nodes2 = Nodes[which(Nodes$phylum_type == "2"), ]
Nodes2 = Nodes2[sort(Nodes2$re_abundance, decreasing = T, index.return = T)$ix, ]
Nodes3 = Nodes[which(Nodes$phylum_type == "3"), ]
Nodes3 = Nodes3[sort(Nodes3$re_abundance, decreasing = T, index.return = T)$ix, ]
Nodes4 = Nodes[which(Nodes$phylum_type == "4"), ]
Nodes4 = Nodes4[sort(Nodes4$re_abundance, decreasing = T, index.return = T)$ix, ]

#----- for order of circle plot -----#
# change 8 and 9, 28 and 32
order = c(30, 31, 28, 33 : 40, 1 : 7, 9, 8, 10 : 27, 32, 29)
Vertices = rbind(Nodes1, Nodes2, Nodes3, Nodes4)
Vertices = Vertices[order, ]
row.names(Vertices) = NULL

#---------------------------------- plot network for CARE -------------------------------------#
Care_lean_net = graph_from_data_frame(d = Care_lean_Edges, vertices = Vertices, directed = FALSE)
Care_obese_net = graph_from_data_frame(d = Care_obese_Edges, vertices = Vertices, directed = FALSE)

color_names = c("#2ca25f", "#e34a33", "#FFFF99", "#BEAED4", "#FDC086", "#9ecae1")

# Run two times for each graph
V(Care_lean_net)$color = (color_names[3 : 6])[as.numeric(V(Care_lean_net)$phylum_type)]
V(Care_lean_net)$size = V(Care_lean_net)$re_abundance + 3
 E(Care_lean_net)$width = 1.5
#E(Care_lean_net)$width = as.numeric(E(Care_lean_net)$weight)
 pdf("care_lean.pdf", width = 9, height = 9)
plot(Care_lean_net, edge.color = (color_names[1 : 2])[(E(Care_lean_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6], 
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()

V(Care_obese_net)$color = (color_names[3 : 6])[as.numeric(V(Care_obese_net)$phylum_type)]
V(Care_obese_net)$size = V(Care_obese_net)$re_abundance + 3
E(Care_obese_net)$width = 1.5
#E(Care_obese_net)$width = as.numeric(E(Care_obese_net)$weight)
pdf("care_obese.pdf", width = 9, height = 9)
plot(Care_obese_net, edge.color = (color_names[1 : 2])[(E(Care_obese_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()


#---------------------------------- plot network for CD-trace -------------------------------------#
CD_lean_net = graph_from_data_frame(d = CD_lean_Edges, vertices = Vertices, directed = FALSE)
CD_obese_net = graph_from_data_frame(d = CD_obese_Edges, vertices = Vertices, directed = FALSE)

color_names = c("#2ca25f", "#e34a33", "#FFFF99", "#BEAED4", "#FDC086", "#9ecae1")

# Run two times for each graph
V(CD_lean_net)$color = (color_names[3 : 6])[as.numeric(V(CD_lean_net)$phylum_type)]
V(CD_lean_net)$size = V(CD_lean_net)$re_abundance + 3
E(CD_lean_net)$width = 1.5
#E(CD_lean_net)$width = as.numeric(E(CD_lean_net)$weight)
pdf("CD_lean.pdf", width = 9, height = 9)
plot(CD_lean_net, edge.color = (color_names[1 : 2])[(E(CD_lean_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()

V(CD_obese_net)$color = (color_names[3 : 6])[as.numeric(V(CD_obese_net)$phylum_type)]
V(CD_obese_net)$size = V(CD_obese_net)$re_abundance + 3
E(CD_obese_net)$width = 1.5
#E(CD_obese_net)$width = as.numeric(E(CD_obese_net)$weight)
pdf("CD_obese.pdf", width = 9, height = 9)
plot(CD_obese_net, edge.color = (color_names[1 : 2])[(E(CD_obese_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()



#---------------------------------- plot network for face -------------------------------------#
face_lean_net = graph_from_data_frame(d = face_lean_Edges, vertices = Vertices, directed = FALSE)
face_obese_net = graph_from_data_frame(d = face_obese_Edges, vertices = Vertices, directed = FALSE)

color_names = c("#4daf4a", "#e41a1c", "#377eb8", "#984ea3", "#a65628", "#ffff33")

# Run two times for each graph
V(face_lean_net)$color = (color_names[3 : 6])[as.numeric(V(face_lean_net)$phylum_type)]
V(face_lean_net)$size = V(face_lean_net)$re_abundance + 3
E(face_lean_net)$width = 1.5
#E(face_lean_net)$width = as.numeric(E(face_lean_net)$weight)
pdf("face_lean.pdf", width = 9, height = 9)
plot(face_lean_net, edge.color = (color_names[1 : 2])[(E(face_lean_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()

V(face_obese_net)$color = (color_names[3 : 6])[as.numeric(V(face_obese_net)$phylum_type)]
V(face_obese_net)$size = V(face_obese_net)$re_abundance + 3
E(face_obese_net)$width = 1.5
#E(face_obese_net)$width = as.numeric(E(face_obese_net)$weight)

pdf("face_obese.pdf", width = 9, height = 9)
plot(face_obese_net, edge.color = (color_names[1 : 2])[(E(face_obese_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()



#---------------------------------- plot network for gCoda -------------------------------------#
gCoda_lean_net = graph_from_data_frame(d = gCoda_lean_Edges, vertices = Vertices, directed = FALSE)
gCoda_obese_net = graph_from_data_frame(d = gCoda_obese_Edges, vertices = Vertices, directed = FALSE)

color_names = c("#2ca25f", "#e34a33", "#FFFF99", "#BEAED4", "#FDC086", "#9ecae1")

# Run two times for each graph
V(gCoda_lean_net)$color = (color_names[3 : 6])[as.numeric(V(gCoda_lean_net)$phylum_type)]
V(gCoda_lean_net)$size = V(gCoda_lean_net)$re_abundance + 3
E(gCoda_lean_net)$width = 1.5
#(gCoda_lean_net)$width = as.numeric(E(gCoda_lean_net)$weight)
pdf("gcoda_lean.pdf", width = 9, height = 9)
plot(gCoda_lean_net, edge.color = (color_names[1 : 2])[(E(gCoda_lean_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()

V(gCoda_obese_net)$color = (color_names[3 : 6])[as.numeric(V(gCoda_obese_net)$phylum_type)]
V(gCoda_obese_net)$size = V(gCoda_obese_net)$re_abundance + 3
E(gCoda_obese_net)$width = 1.5
#E(gCoda_obese_net)$width = as.numeric(E(gCoda_obese_net)$weight)
pdf("gcoda_obese.pdf", width = 9, height = 9)
plot(gCoda_obese_net, edge.color = (color_names[1 : 2])[(E(gCoda_obese_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()





#---------------------------------- plot network for SPIEC-EASI -------------------------------------#
SE_gl_lean_net = graph_from_data_frame(d = SE_gl_lean_Edges, vertices = Vertices, directed = FALSE)
SE_gl_obese_net = graph_from_data_frame(d = SE_gl_obese_Edges, vertices = Vertices, directed = FALSE)

color_names = c("#2ca25f", "#e34a33", "#FFFF99", "#BEAED4", "#FDC086", "#9ecae1")

# Run two times for each graph
V(SE_gl_lean_net)$color = (color_names[3 : 6])[as.numeric(V(SE_gl_lean_net)$phylum_type)]
V(SE_gl_lean_net)$size = V(SE_gl_lean_net)$re_abundance + 3
E(SE_gl_lean_net)$width = 1.5
#E(SE_gl_lean_net)$width = as.numeric(E(SE_gl_lean_net)$weight)
pdf("SE_lean.pdf", width = 9, height = 9)
plot(SE_gl_lean_net, edge.color = (color_names[1 : 2])[(E(SE_gl_lean_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()

V(SE_gl_obese_net)$color = (color_names[3 : 6])[as.numeric(V(SE_gl_obese_net)$phylum_type)]
V(SE_gl_obese_net)$size = V(SE_gl_obese_net)$re_abundance + 3
E(SE_gl_obese_net)$width = 1.5
#E(SE_gl_obese_net)$width = as.numeric(E(SE_gl_obese_net)$weight)
pdf("SE_obese.pdf", width = 9, height = 9)
plot(SE_gl_obese_net, edge.color = (color_names[1 : 2])[(E(SE_gl_obese_net)$sign == "-1") + 1],
     layout = layout_in_circle, edge.curved = .40, vertex.label.family= "sans", 
     vertex.label.cex=.75, vertex.label.color = "black", vertex.frame.color = NA, margin = c(-0.3, 0, -0.3, 0.3))
legend(x = 1.1, y = -0.2, c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria"), pch = 21, pt.bg = color_names[3 : 6],
       pt.lwd = 0, bty = "n", pt.cex = 1.0, cex = 0.75, ncol = 1, x.intersp = 0.5, y.intersp = 0.7)
dev.off()




