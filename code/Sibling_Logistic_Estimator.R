logit <- function(x){
  return( exp(x)/(1 + exp(x)))
}

library(Matrix)


format_data <- function(pheno, F_ind, Z = NULL){
  Z_tilde = NULL
  
  F_ind = match(F_ind, unique(F_ind))
  # Number of families
  K = max(F_ind)
  # Number of individuals
  N = length(F_ind)
  
  size_dic <- as.integer(table(F_ind))
  F_size <- size_dic[F_ind]
  
  family2ind = list()
  for(i in 1:K){
    family2ind[[i]] = which(F_ind == i)
  }
  
  if(1 %in% table(F_ind) ){
    warning("Some families only have 1 individual")
  }
  pheno_tilde <- pheno
  family_means <- tapply(pheno_tilde, F_ind, mean)
  pheno_tilde <- pheno_tilde - family_means[F_ind]
  
  if(! is.null(Z)){
    Z = as.matrix(model.matrix( ~ -1 + ., data = Z))
    Z_tilde = Z
    for(k in 1:ncol(Z)){
      kth_mean = tapply(Z[, k], F_ind, mean)
      Z_tilde[, k] = Z[, k] - kth_mean[F_ind]
    }
  }
  
  dim_Z = 0
  if(! is.null(Z)){
    dim_Z = length(Z) / N
  } 
  
  return(list(Y = pheno, Y_tilde = as.vector(pheno_tilde), F_ind = F_ind, Z = Z, Z_tilde = Z_tilde,
              size_dic = size_dic, F_size = F_size, K = max(F_ind),  N = length(F_ind), dim_Z = dim_Z,
              family2ind = family2ind))
}


calibrated_logistic_estimator <- function(X, data, alpha_ext, alpha_ext_var, N_ext,
                                 overlap_ratio = 0){
  with(data, {
    # Adjust external variance
    C11 = (alpha_ext_var / N_ext) * N

    # Process data    
    X = matrix(X, N, 1)
    family_effect = tapply(X, F_ind, sum)
    f = family_effect[F_ind]
    XZ = cbind(X, Z)
    
    # Compute mis-specified model
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
    
    # Efficiently create diagonal matrices using Matrix package
    coef_cor <- Diagonal(x = predict_cor * (1 - predict_cor))
    coef_mis <- Diagonal(x = predict_mis * (1 - predict_mis))
    
    # Efficient matrix multiplications and scaling
    D2 <- crossprod(XZ_tilde, coef_cor %*% XZ_tilde) / N
    D3 <- crossprod(XZ, coef_mis %*% XZ) / N
    


    # Compute the C_{22}, C_{33}, C_{23}

    M = max(F_size)
    
    correct_var = list()
    incorrect_var = list()
    correct_incorrect_cov = list()
    
    dim_XZ = dim(XZ)[2]
    # Consider each family size separately
    if(M == 2){
      prod = (Y - predict_mis) * XZ
      prod_tilde =  (Y - predict_cor) * XZ_tilde
      
      
      num_F_ind_levels <- length(unique(F_ind))
      est <- matrix(0, nrow = num_F_ind_levels, ncol = (dim_XZ + 1) + dim_XZ)
      
      # Compute the first part of est
      for(i in 1:(dim_XZ + 1)){
        est[, i] <- tapply(prod_tilde[, i], F_ind, sum)
      }
      
      # Compute the second part of est
      for(i in 1:dim_XZ){
        est[, (dim_XZ + 1 + i)] <- tapply(prod[, i], F_ind, sum)
      }
      
      # Calculate the covariance of est
      temp <- cov(est)
      
      
      correct_var = temp[1:(dim_XZ + 1), 1:(dim_XZ + 1)]
      
      # This is for C33
      incorrect_var = temp[(dim_XZ+2): (2 * dim_XZ + 1), (dim_XZ+2): (2 * dim_XZ + 1)]
      
      # This is for C23
      correct_incorrect_cov = temp[1:(dim_XZ + 1), (dim_XZ+2):(2 * dim_XZ + 1)]
      
      C22 =  correct_var * K / N
      C33 = incorrect_var * K / N
      C23 =  correct_incorrect_cov  * K / N
      
      
    } else{
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
      C22 =  Reduce("+", correct_var[size_dic[1:K]]) / N
      C33 = Reduce("+", incorrect_var[size_dic[1:K]]) / N
      C23 = Reduce("+", correct_incorrect_cov[size_dic[1:K]]) / N
    }




    
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
    
    sub_ind = sample_indices(F_ind, max(F_ind))
    
    int_mis = glm(Y[sub_ind] ~ XZ[sub_ind, ], family = binomial(link = "logit"))
    alpha_int_single = summary(int_mis)$coefficients[2,1]
    alpha_int_var_single = summary(int_mis)$coefficients[2,2]^2
    
    return(list(beta_cal  = beta_cal, beta_cal_var = (beta_cal_var), beta_int = beta_int, 
                beta_int_var = (C22/N), alpha_int = alpha_int, alpha_int_var = (C33)/N, 
                alpha_int_single = alpha_int_single,alpha_int_var_single = alpha_int_var_single))
   } )
}

sample_indices <- function(vec, K) {
  sampled_indices <- sapply(1:K, function(k) {
    indices <- which(vec == k)
    if (length(indices) > 0) {
      sample(indices, 1)
    } else {
      NA  # Handle the case where there are no elements equal to k
    }
  })
  return(sampled_indices)
}

