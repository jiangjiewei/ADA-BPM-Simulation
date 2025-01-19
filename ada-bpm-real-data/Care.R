###################### Composition Adaptive Regularized Estimation (CARE) ######################
#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
# Care_col: 
#     computes the estimator of omega_j for the jth column under different tuning parameters.
# CV_Care_col: 
#     returns the selected tuning parameter by cross-validation for the jth column.
# Care_est: 
#     computes the estimator of Omega under the selected tuning parameter.
# Care_sp: 
#     computes the estimators of Omega under different tuning parameters.
#----------------------------------------------------------------------------------------------#

Care_col <- function(col_num, X, nlambda, lambda_min){
  #----------------------------------------------------------------------#
  # Input: 
  #       col_num, numeric,
  #           the index of each column; 
  #       X, (n * p) matrix, 
  #           compositional data;
  #       nlambda, numeric, 
  #           the number of tuning parameters;
  #       lambda_min, numeric, 
  #           the minimal tuning parameter.
  # Output:
  #       omega_hat_mat, (p * nlambda) matrix, 
  #           the estimator of omega_j under different tuning parameters;
  #       lambda_vec, (nlambda * 1) vector,
  #           the sequence of tuning parameters.
  #----------------------------------------------------------------------#
  
  n <- nrow(X)
  p <- ncol(X)
  G <- diag(p) - matrix(1, p, p) / p
  Z <- log(X) %*% G
  Sigma_c_hat <- cov(Z) * (1 - 1 / n) 
  A <- cbind(cbind(rbind(Sigma_c_hat, -Sigma_c_hat), rbind(-Sigma_c_hat, Sigma_c_hat)), 
             diag(rep(1, 2 * p)))
  b <- c(diag(p)[, col_num] - rep(1, p) / p, -diag(p)[, col_num] + rep(1, p) / p)
  c <- c(rep(-1, 2 * p), rep(0, 2 * p))
  c_bar <- rep(0, 4 * p)
  b_bar <- rep(1, 2 * p)
  B_init <- seq(2 * p, 4 * p - 1)
  output <- PSM_solver(A, b, b_bar, c, c_bar, B_init, max_it = nlambda, 
                       lambda_threshold = lambda_min)
  lambda_vec <- output$lambda
  solution <- as.matrix(output$beta)[1 : (2 * p), ]
  omega_hat_mat<- solution[1 : p, ] - solution[(p + 1) : (2 * p), ]
  return(list(omega_hat_mat = omega_hat_mat, lambda_vec = lambda_vec))
}


CV_Care_col = function(col_num, X, lambda_min, lambda_vec, ratio, nsplit){
  #------------------------------------------------------------#
  # Input: 
  #       col_num, numeric,
  #           the index of each column; 
  #       X, (n * p) matrix, 
  #           compositional data;
  #       lambda_min, numeric, 
  #           the minimal tuning parameter;
  #       lambda_vec, (nlambda * 1) vector, 
  #           the sequence of tuning parameters;
  #       ratio, numeric,
  #           the ratio of training sample in the whole sample;
  #       nsplit, numeric, 
  #           the number of sample splitting.
  # Output:
  #       which.min(error), numeric, 
  #           the index of the selected tuning parameter.
  #------------------------------------------------------------#
  
  n <- nrow(X)
  p <- ncol(X)
  G <- diag(p) - matrix(1, p, p) / p
  nlambda <- length(lambda_vec)
  error <- rep(0, nlambda)
  for(i in 1 : nsplit){
    train <- sort(sample(1 : n, round(ratio * n)))
    test <- setdiff(1 : n, train)
    output <- Care_col(col_num, X[train, ], nlambda * 2, lambda_min) 
    index <- apply(abs(outer(lambda_vec, output$lambda_vec, "-")), 1, which.min)
    omega_hat_mat_train <- (output$omega_hat_mat)[, index]
    Sigma_c_hat_test <- cov(log(X[test, ]) %*% G) * (1 - 1 / length(test))
    Loss <- diag(t(omega_hat_mat_train) %*% Sigma_c_hat_test %*% omega_hat_mat_train) / 2 - 
            t(omega_hat_mat_train) %*% (diag(p)[, col_num] - rep(1, p) / p)
    error <- error + as.numeric(Loss) / nsplit
  }
  return(which.min(error))
}


Care_est <- function(X, nlambda, lambda_min, ratio, nsplit){
  #----------------------------------------------------------------------#
  # Input: 
  #       X, (n * p) matrix, 
  #           compositional data;
  #       nlambda, numeric, 
  #           the number of tuning parameters;
  #       lambda_min, numeric, 
  #           the minimal tuning parameter;
  #       ratio, numeric,
  #           the ratio of training sample in the whole sample;
  #       nsplit, numeric, 
  #           the number of sample splitting.
  # Output:
  #       Omega_hat, (p * p) matrix, 
  #           the estimator of Omega under the selected tuning parameter;
  #       lambda_op, (p * 1) vector,
  #           the selected tuning parameter for each column.
  #----------------------------------------------------------------------#
  
  p <- ncol(X)
  Omega_hat <- matrix(0, p, p)
  lambda_op <- rep(0, p)
  for(j in 1 : p){
    output <- Care_col(col_num = j, X, nlambda, lambda_min)
    opt <- CV_Care_col(col_num = j, X, lambda_min, lambda_vec = output$lambda_vec, ratio, nsplit)
    Omega_hat[, j] <- (output$omega_hat_mat)[, opt]
    lambda_op[j] <- (output$lambda_vec)[opt]
  }
  Omega_hat <- Omega_hat * (abs(Omega_hat) <= abs(t(Omega_hat))) + 
               t(Omega_hat) * (abs(Omega_hat) > abs(t(Omega_hat))) 
  return(list(Omega_hat = Omega_hat, lambda_op = lambda_op))
}


Care_sp <- function(X, nlambda, lambda_min){
  #----------------------------------------------------------------------#
  # Input: 
  #       X, (n * p) matrix, 
  #           compositional data;
  #       nlambda, numeric, 
  #           the number of tuning parameters;
  #       lambda_min, numeric, 
  #           the minimal tuning parameter.
  # Output:
  #       path, (nlambda * 1) list, 
  #           the estimators of Omega under different tuning parameters.
  #----------------------------------------------------------------------#
  
  p <- ncol(X)
  omega_hat_mat <- NULL
  for(j in 1 : p){
    output <- Care_col(col_num = j, X, nlambda, lambda_min)
    omega_hat_mat <- rbind(omega_hat_mat, output$omega_hat_mat)
  }
  path <- list()
  for(i in 1 : nlambda){
    Omega_hat <- matrix(omega_hat_mat[, i], nrow = p, ncol = p)
    path[[i]] <- Omega_hat * (abs(Omega_hat) <= abs(t(Omega_hat))) + 
                 t(Omega_hat) * (abs(Omega_hat) > abs(t(Omega_hat))) 
  }
  return(path)
}

