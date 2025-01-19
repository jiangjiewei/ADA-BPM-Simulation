#----------------------- Construct network dataset -----------------------#
Edges_data = function(Omega_0_hat, Nodes){
  #-----------------------------------#
  # Input: 
  #       Omega_0_hat, (p * p) matrix, 
  #           the estimator of Omega;
  #       Nodes, data frame, 
  #           the names of nodes. 
  #
  # Output:
  #       Edges, data frame, 
  #            network dataset.
  #-----------------------------------#
  
  p = ncol(Omega_0_hat)
  from = NULL; to = NULL; sign = NULL; strength = NULL
  genus = (as.matrix(Nodes))[, 1]
  
  for(i in 1 : (p - 1)){
    candidate = (i + 1) : p
    temp = which(Omega_0_hat[candidate, i] != 0)
    strength = c(strength, round((Omega_0_hat[candidate, i])[temp], 5))
    
    # edge names (node to node)
    from = c(from, rep(genus[i], length(temp)))
    to = c(to, genus[candidate[temp]])
    
    # whether edges between connected nodes are positive or negative 
    # 1 positive   v.s.  -1 negative 
    sign = c(sign, ifelse((Omega_0_hat[candidate, i])[temp] > 0, 1, -1))
    
    # transform the strength for plot
    # weight = round(log((abs(strength) * 1000))) + 3
    weight = round((abs(strength) * 100)) + 2
    weight[which(weight >= 20)] = 20
    
  }
  Edges = data.frame(cbind(from, to, sign, strength, weight))
  return(Edges)
}


#--------------------- Compute positive or negative edges ---------------------#
Strength = function(Omega_0_hat){
  #---------------------------------------#
  # Input: 
  #       Omega_0_hat, (p * p) matrix, 
  #           the estimator of Omega.
  # Output:
  #       Num_positive, numeric,
  #           the number of positive edges;
  #       Num_negative, numeric,
  #           the number of negative edges.
  #---------------------------------------#
  
  Lower_Omega_0_hat = Omega_0_hat[lower.tri(Omega_0_hat)]
  Num_positive = length(which(Lower_Omega_0_hat > 0))
  Num_negative = length(which(Lower_Omega_0_hat < 0))
  return(c(Num_positive, Num_negative))
}


#------------------ Compute network stability and edge replicates ------------------#
Subsampling = function(Omega_0_hat_list, compdata, pro, num){
  #----------------------------------------------------------------------#
  # Input: 
  #       Omega_0_hat_list, list, 
  #           the estimators of Omega obtained by different methods;
  #       compdata, data matrix,
  #           the compositional dara;
  #       pro, numeric,
  #           the ratio of subsampling;
  #       num, numeric,
  #           the replicates of subsampling.
  # Output:
  #       Stability, 5 * 1 vector,
  #           the network stability ;
  #       Care_replicate, CD_replicate, gCoda_replicate, SE_gl_replicate, 
  #           edge replicates for different methods.
  #-----------------------------------------------------------------------#
  
  
  Care_Omega_0_hat = Omega_0_hat_list[[1]]
  CD_Omega_0_hat = Omega_0_hat_list[[2]]
  face_Omega_0_hat=Omega_0_hat_list[[3]]
  gCoda_Omega_0_hat = Omega_0_hat_list[[4]] 
  SE_gl_Omega_0_hat = Omega_0_hat_list[[5]]
  
  Care_lower = Care_Omega_0_hat[lower.tri(Care_Omega_0_hat)]
  CD_lower = CD_Omega_0_hat[lower.tri(CD_Omega_0_hat)]
  face_lower = face_Omega_0_hat[lower.tri(face_Omega_0_hat)]
  gCoda_lower = gCoda_Omega_0_hat[lower.tri(gCoda_Omega_0_hat)]
  SE_gl_lower = SE_gl_Omega_0_hat[lower.tri(SE_gl_Omega_0_hat)]
  
  Care_nedge = length(which(Care_lower != 0)); Care_replicate = rep(0, Care_nedge)
  CD_nedge = length(which(CD_lower != 0)); CD_replicate = rep(0, CD_nedge)
  face_nedge = length(which(face_lower != 0)); face_replicate = rep(0, face_nedge)
  gCoda_nedge = length(which(gCoda_lower != 0)); gCoda_replicate = rep(0, gCoda_nedge)
  SE_gl_nedge = length(which(SE_gl_lower != 0)); SE_gl_replicate = rep(0, SE_gl_nedge)

  n = nrow(compdata)
  p = ncol(compdata)
  G = diag(p) - matrix(1, p, p) / p
  
  Care_stab = rep(0, num)
  CD_stab = rep(0, num)
  face_stab = rep(0, num)
  gCoda_stab = rep(0, num)
  SE_gl_stab = rep(0, num)
  
  for(i in 1 : num){
    set.seed(i)
    sample_id = sort(sample(1 : n, round(n * pro)))
    new_compdata = compdata[sample_id, ]
    new_compdata_clr = log(new_compdata) %*% G
    Nlambda = 50
    
    FsigmaF<-cov(new_compdata_clr)*(1-1/n);
    lam.max <- max(max(FsigmaF - diag(p)), -min(FsigmaF - diag(p)));
    lambda.min.ratio4=0.011 #0.01
    lam.min4 <- lambda.min.ratio4 * lam.max;
    lam4 <- exp(seq(log(lam.max), log(lam.min4), length = Nlambda));
    
    #------------------------- CARE -------------------------#
    Care_sub_output = Care_est(new_compdata, nlambda = Nlambda, lambda_min = 0.01, ratio = 0.75, nsplit = 10)
    Care_sub_Omega_0_hat = Care_sub_output$Omega_hat
    Care_sub_lower = Care_sub_Omega_0_hat[lower.tri(Care_sub_Omega_0_hat)]
    Care_temp = ifelse(Care_lower != 0, 1, 0) * ifelse(Care_sub_lower != 0, 1, 0)
    Care_replicate = Care_temp[which(Care_lower != 0)] + Care_replicate
    Care_stab[i] = sum(Care_sub_lower[Care_lower != 0] != 0) / Care_nedge
    
    #----------------------- CD-trace -------------------------#
    CD_output = CompDtrace(new_compdata_clr, nlambda = Nlambda, lambda.min.ratio = 0.01)
    CD_sub_Omega_0_hat = (CD_output$Theta)[[which.min(CD_output$bic)]]
    CD_sub_lower = CD_sub_Omega_0_hat[lower.tri(CD_sub_Omega_0_hat)]
    CD_temp = ifelse(CD_lower != 0, 1, 0) * ifelse(CD_sub_lower != 0, 1, 0)
    CD_replicate = CD_temp[which(CD_lower != 0)] + CD_replicate
    CD_stab[i] = sum(CD_sub_lower[CD_lower != 0] != 0) / CD_nedge
    
    #--------------------------- face ------------------------------#
    face_output=ADA_BPM(new_compdata_clr,lambda =lam4, gamma=0.999)
    face_sub_Omega_0_hat=face_output$Omega[[which.min(face_output$bic)]]
    face_sub_lower = face_sub_Omega_0_hat[lower.tri(face_sub_Omega_0_hat)]
    face_temp = ifelse(face_lower != 0, 1, 0) * ifelse(face_sub_lower != 0, 1, 0)
    face_replicate = face_temp[which(face_lower != 0)] + face_replicate
    face_stab[i] = sum(face_sub_lower[face_lower != 0] != 0) / face_nedge
    
    #------------------------- gCoda --------------------------#
    gCoda_sub_Omega_0_hat = gcoda(new_compdata_clr, nlambda = Nlambda, lambda.min.ratio = 0.01)$opt.icov
    gCoda_sub_lower = gCoda_sub_Omega_0_hat[lower.tri(gCoda_sub_Omega_0_hat)]
    gCoda_temp = ifelse(gCoda_lower != 0, 1, 0) * ifelse(gCoda_sub_lower != 0, 1, 0)
    gCoda_replicate = gCoda_temp[which(gCoda_lower != 0)] + gCoda_replicate
    gCoda_stab[i] = sum(gCoda_sub_lower[gCoda_lower != 0] != 0) / gCoda_nedge
    
    #----------------------- S-E(glasso) --------------------------#
    SE_gl_output = sparseiCov(new_compdata_clr, method = "glasso", lambda.min.ratio = 0.01, nlambda = Nlambda)
    SE_gl_sub_Omega_0_hat = as.matrix(huge.select(SE_gl_output, criterion = "stars", verbose = FALSE)$opt.icov)
    SE_gl_sub_lower = SE_gl_sub_Omega_0_hat[lower.tri(SE_gl_sub_Omega_0_hat)]
    SE_gl_temp = ifelse(SE_gl_lower != 0, 1, 0) * ifelse(SE_gl_sub_lower != 0, 1, 0)
    SE_gl_replicate = SE_gl_temp[which(SE_gl_lower != 0)] + SE_gl_replicate
    SE_gl_stab[i] = sum(SE_gl_sub_lower[SE_gl_lower != 0] != 0) / SE_gl_nedge
  }
  Stability = c(mean(Care_stab), mean(CD_stab), mean(face_stab), mean(gCoda_stab), mean(SE_gl_stab))
  return(list(Stability, Care_replicate, CD_replicate, face_replicate, gCoda_replicate, SE_gl_replicate))
}


