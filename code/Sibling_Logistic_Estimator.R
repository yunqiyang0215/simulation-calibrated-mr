##' Compute the calibrated estimator
##' @param Y: a vector contains response variable values
##' @param X: a vector contains transmitted allele values
##' @param F_ind: family index (from 1 to K)
##' @param alpha_ext: estimate of alpha from external data
##' @param alpha_ext_var: variance of alpha from external data
##' @param N_ext: number of samples in external data
##' @param overlap_ratio: proportion of internal data from external data
##' 
calibrated_logistic_estimator <- function(Y, X, F_ind, alpha_ext, alpha_ext_var, N_ext, overlap_ratio = 0){
  
  # Number of families
  K = max(F_ind)
  # Number of individuals
  N = length(X)
  
  # Compute family size for each individual
  size_dic = rep(0, K)
  F_size = rep(0, N)
  for(i in 1:N){
    size_dic[F_ind[i]] = size_dic[F_ind[i]]  + 1
  }
  
  for(i in 1:N){
    F_size[i] = size_dic[F_ind[i]]
  }
  
  family_effect = rep(0, K)
  for(i in 1:N){
    family_effect[F_ind[i]] = family_effect[F_ind[i]] + X[i]
  }
  
  f = family_effect[F_ind]
  
  # Adjust external variance
  C11 = (alpha_ext_var / N_ext) * N
  
  # Compute the incorrect model for internal data
  
  int_mis = summary(glm(Y ~ X, family = binomial(link = "logit")))
  
  alpha_int = int_mis$coefficients[2]
  predict_mis = logit( int_mis$coefficients[1] +  int_mis$coefficients[2] * X)
  # Compute the correct model for internal data
  int_cor = summary(glm(Y ~ X + f, family = binomial(link = "logit")))
  beta_int = int_cor$coefficients[2,1]
  predict_cor = logit( int_cor$coefficients[1] +  int_cor$coefficients[2] * X 
                       + int_cor$coefficients[3] * f)
  
  # Compute the variance
  
  X = cbind(1, X)
  X_tilde =cbind(X, f)
  
  D2 =  matrix(rep(0, 9), 3, 3)
  D3 = matrix(rep(0, 4), 2, 2)
  
  for(i in 1:N){
    D2 = D2 + predict_cor[i] * (1 - predict_cor[i]) *( X_tilde[i, ] %*% t(X_tilde[i, ]))
    D3 = D3 + predict_mis[i] * (1 - predict_mis[i]) *( X[i, ] %*% t(X[i, ]) )
  }
  
  D2 = D2 /  N
  D3 = D3 / N
  
  # Compute the C_{22}, C_{33}, C_{23}
  
  M = max(F_size)
  
  correct_var = list()
  incorrect_var = list()
  correct_incorrect_cov = list()
  
  # Consider each family size separately
  for(m in 2:M){
    
    # Only look at family of size m
    X_sub = X[F_size == m, ]
    predict_mis_sub = predict_mis[F_size == m]
    
    X_tilde_sub = X_tilde[F_size == m, ]
    predict_cor_sub = predict_cor[F_size == m]
    
    Y_sub = Y[F_size == m]
    
    prod = (Y_sub - predict_mis_sub) * X_sub 
    prod_tilde =  (Y_sub - predict_cor_sub) * X_tilde_sub
    
    
    temp1 = (tapply(prod_tilde[,1], F_ind[F_size == m], sum))
    temp2 = (tapply(prod_tilde[,2], F_ind[F_size == m], sum))
    temp3 = (tapply(prod_tilde[,3], F_ind[F_size == m], sum))
    temp4 = tapply(prod[,1], F_ind[F_size == m], sum)
    temp5 = tapply(prod[,2], F_ind[F_size == m], sum)
    
    temp = cov(cbind(temp1, temp2, temp3, temp4, temp5))
    
    # This is for C22
    
    correct_var[[m]] = temp[1:3, 1:3]
    # This is for C33
    
    incorrect_var[[m]] = temp[4:5, 4:5]
    # This is for C23
    correct_incorrect_cov[[m]] = temp[1:3, 4:5]
  }  
  
  
  C22 = 0
  C33 = 0
  C23 = 0
  
  for(i in 1:K){
    C22 = C22 + correct_var[[F_size[i]]]
    C33 = C33 + incorrect_var[[F_size[i]]]
    C23 = C23 + correct_incorrect_cov[[F_size[i]]]
  }
  
  
  C22 = ( solve(D2) %*% C22 %*% solve(D2) )[2,2] / N
  C33 = (solve(D3) %*% C33 %*% solve(D3) )[2,2] / N
  C23 = (solve(D2) %*% C23 %*% solve(D3))[2,2] / N
  
  # Compute C12, C13
  
  C12 = overlap_ratio * C23 * N / N_ext
  C13 = overlap_ratio * C33 * N / N_ext
  
  # Covariance matrix of (alpha_ext - alpha_int, beta_int - beta)
  result_cov = rbind(c(C11 - 2 * C13+ C33, C12 - C23), c(C12 - C23, C22))
  # Compute the final calibrated estimator
  beta_cal = beta_int +  (C23 - C12) / ( C11 + C33 - 2 * C13) * (alpha_ext - alpha_int)
  beta_cal_var = (C22 - (C23 - C12)^2 / (C11 + C33 - 2 * C13) ) / N
  return(list(beta_cal  = beta_cal, beta_cal_var = (beta_cal_var), beta_int = beta_int, 
              beta_int_var = (C22/N)))
}



logit <- function(x){
  return( exp(x)/(1 + exp(x)))
}
