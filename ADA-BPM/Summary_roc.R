########################### Summary of Support Recovery ############################
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
# ROC_PR: 
#     computes the true positive and false positive rates of an estimator of Omega.
# summary_ROC_PR: 
#     computes the true positive and false positive rates for the solution path.
#----------------------------------------------------------------------------------#

ROC_PR <- function(pred, true){
  #-----------------------------------#
  # Input: 
  #       pred, (p * p) matrix, 
  #           the estimator of Omega;
  #       true, (p * p) matrix, 
  #           true Omega. 
  #
  # Output:
  #       TPR, numeric, 
  #            true positive rate;
  #       FPR, numeric, 
  #            false positive rate.
  #-----------------------------------#
  
  Upper_true <- true[upper.tri(true)]
  Upper_pred <- pred[upper.tri(pred)]
  TP <- sum(Upper_pred[Upper_true != 0] != 0)
  FP <- sum(Upper_pred[Upper_true == 0] != 0)
  FN <- sum(Upper_pred[Upper_true != 0] == 0)
  TN <- sum(Upper_pred[Upper_true == 0] == 0)
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  return(c(TPR, FPR))
}


summary_ROC_PR = function(path, true){
  #---------------------------------------------------------------------#
  # Input: 
  #       path, (nlambda * 1) list, 
  #           the estimators of Omega under different tuning parameters;
  #       true, (p * p) matrix, 
  #           true Omega. 
  #
  # Output:
  #       result, (nlambda * 2) matrix, 
  #           TPR and FPR under different tuning parameters.
  #---------------------------------------------------------------------# 
  
  result <- data.frame(t(sapply(path, function(pred)ROC_PR(pred, true))))
  colnames(result) <- c('TPR', 'FPR')
  return(result)
}
