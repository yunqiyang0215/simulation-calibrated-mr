calibrated_logistic_estimator <- function(Y, X, F_ind, alpha_ext, alpha_ext_var, N_ext, overlap_ratio = 0,
                                          Z = NULL){
  
  
  # Reindex F_ind
  
  F_ind = match(F_ind, unique(F_ind))
  
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
  if(is.null(Z)){
    dim_Z = 0
  } else{
    dim_Z = dim(Z)[2]
  }
  XZ = cbind(X, Z)
  XZ = matrix(XZ, ncol = 1 + dim_Z)
  
  
  int_mis = glm(Y ~ XZ, family = binomial(link = "logit"))
  alpha_int = as.numeric(int_mis$coefficients[2])
  
  predict_mis = logit(predict(int_mis))
  
  # Compute the correct model for internal data
  int_cor = glm(Y ~ XZ + f, family = binomial(link = "logit"))
  beta_int = as.numeric(int_cor$coefficients[2])
  predict_cor = logit(predict(int_cor))
  
  # Compute the variance
  
  XZ = cbind(1, XZ)
  XZ_tilde =cbind(XZ, f)
  
  coef_cor = diag(predict_cor * (1 - predict_cor))
  coef_mis = diag(predict_mis * (1 - predict_mis))
  
  D2 = t(XZ_tilde) %*% coef_cor %*% XZ_tilde
  D3 = t(XZ) %*% coef_mis %*% XZ
  

  D2 = D2 /  N
  D3 = D3 / N
  
  # Compute the C_{22}, C_{33}, C_{23}
  
  M = max(F_size)
  
  correct_var = list()
  incorrect_var = list()
  correct_incorrect_cov = list()
  
  dim_XZ = dim(XZ)[2]
  # Consider each family size separately
  for(m in 2:M){
    
    # Only look at family of size m
    XZ_sub = XZ[F_size == m, ]
    predict_mis_sub = predict_mis[F_size == m]
    
    XZ_tilde_sub = XZ_tilde[F_size == m, ]
    predict_cor_sub = predict_cor[F_size == m]
    
    Y_sub = Y[F_size == m]
    
    prod = (Y_sub - predict_mis_sub) * XZ_sub 
    prod_tilde =  (Y_sub - predict_cor_sub) * XZ_tilde_sub
    

    # This is for C22
    est = NULL
    
    for(i in 1: (dim_XZ + 1)){
      est = cbind(est, tapply(prod_tilde[,i], F_ind[F_size == m], sum))
    }
    for(i in 1:dim_XZ){
      est = cbind(est,  tapply(prod[,i], F_ind[F_size == m], sum))
    }
    
    temp = cov(est)
    
    
    correct_var[[m]] = temp[1:(dim_XZ + 1), 1:(dim_XZ + 1)]
    
    # This is for C33
    incorrect_var[[m]] = temp[(dim_XZ+2): (2 * dim_XZ + 1), (dim_XZ+2): (2 * dim_XZ + 1)]
    
    # This is for C23
    correct_incorrect_cov[[m]] = temp[1:(dim_XZ + 1), (dim_XZ+2):(2 * dim_XZ + 1)]
  }  
  
  
  C22 = 0
  C33 = 0
  C23 = 0
  
  for(i in 1:K){
    C22 = C22 + correct_var[[F_size[i]]]
    C33 = C33 + incorrect_var[[F_size[i]]]
    C23 = C23 + correct_incorrect_cov[[F_size[i]]]
  }
  
  C22 = C22 / N
  C33 = C33 / N
  C23 = C23 /N
  
  C22 = ( solve(D2) %*% C22 %*% solve(D2) )[2,2] 
  C33 = (solve(D3) %*% C33 %*% solve(D3) )[2,2] 
  C23 = (solve(D2) %*% C23 %*% solve(D3))[2,2] 
  
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
