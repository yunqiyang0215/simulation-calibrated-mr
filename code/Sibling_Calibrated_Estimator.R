library(mvtnorm)


##' Generate (Y, X, F_ind) such that each family has exactly two individuals
##'
##' @param sigma_esp: standard deviation of noise epsilon
##' @param eps_cor: correlation of noise within a family
##' @param sigma_x: standard deviation of X
##' @param sigma_f: standard deviation of f
##' @param rho: correlation between X and f
##' @param beta: true linear model coefficient
##' @param K: number of families
##'
generate_data <- function(sigma_eps = .5, eps_cor = .1, sigma_x =1,
                          sigma_f = 1, rho = .3, beta = 1, K = 1000, seed = 123){
  set.seed(seed)
  # Compute the covariance matrix for (X_{k1}, X_{k2}, f_k)
  Sigma = rbind(c(sigma_x^2, 0,  rho * sigma_x * sigma_f),
                c(0, sigma_x^2,  rho * sigma_x * sigma_f),
                c(sigma_x * sigma_f * rho,  sigma_x * sigma_f * rho, sigma_f^2))

  # Generate X, f
  x1_x2_f = rmvnorm(K, sigma = Sigma)
  x_1 = x1_x2_f[, 1]
  x_2 = x1_x2_f[, 2]
  f = x1_x2_f[, 3]

  # Generate noise
  noise = rmvnorm(K, sigma = rbind( c(sigma_eps^2, sigma_eps^2 * eps_cor), c(sigma_eps^2 * eps_cor, sigma_eps^2) ))

  # Generate Y
  y_1 = beta * x_1 + f + noise[, 1]
  y_2 = beta * x_2 + f + noise[, 2]
  Y = c(y_1, y_2)
  X = c(x_1, x_2)

  F_ind = c(1:K, 1:K)

  return(list(Y = Y, X = X, F_ind = F_ind))

}

##‘ Sample one sibling from each family
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

library(dplyr)


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

##' Compute the calibrated estimator
##' @param X: a vector that contains transmitted allele values
##' @param data: output from format_data
##' @param alpha_ext: estimate of alpha from external data
##' @param alpha_ext_var: variance of alpha from external data
##' @param N_ext: number of samples in external data
##' @param overlap_ratio: proportion of internal data from external data
##'
calibrated_estimator <- function(X, data, alpha_ext, alpha_ext_var, N_ext,
                                 overlap_ratio = 0){
  with(data, {
    # Adjust external variance
    C11 = (alpha_ext_var / N_ext) * N
    
    # Compute correct model
    
    X = matrix(X, N, 1)
    X_tilde = X - matrix(tapply(X, F_ind, mean)[F_ind], N, 1)
    
    XZ = cbind(X, Z)
    XZ_tilde = cbind(X_tilde, Z_tilde)
    
    
    correct_int = lm(Y_tilde ~ -1 +  XZ_tilde)
    beta_int = summary(correct_int)$coefficients[1, 1]
    
    # Compute the incorrect model for internal data
    
    incorrect_int = lm( Y ~  XZ)
    alpha_int = summary(incorrect_int)$coefficients[2, 1]
    # Compute the C_{22}, C_{33}, C_{23}
    
    M = max(F_size)
    Cs = list()
    # Consider each family size separately
    if(M == 2){
      resid = resid(summary(incorrect_int))
      
      resid_tilde = resid(summary(correct_int))
      
      
      C = cbind(XZ_tilde * resid_tilde, resid, XZ * resid)
      
      C = t(sapply(1:K,
                   function(x) colSums(C[family2ind[[x]], ])))
      
      
      Cs[[2]] = cov(C)
      
    }
    else {
    for(m in 2:M){
      
      # Only look at family of size m
      XZ_sub = XZ[F_size == m, ]
      resid_sub = resid(summary(incorrect_int))[F_size == m]
      
      XZ_tilde_sub = XZ_tilde[F_size == m, ]
      resid_tilde_sub = resid(summary(correct_int))[F_size == m]
      
      
      ind_sub = F_ind[F_size == m]
      
      C = cbind(XZ_tilde_sub * resid_tilde_sub, resid_sub, XZ_sub * resid_sub)
      
      C = t(sapply(unique(ind_sub),
                   function(x) colSums(C[family2ind[[x]], ])))
      
      
      Cs[[m]] = cov(C)
      
    }}
    

    final_C <- Reduce("+", Cs[size_dic[1:K]])
    

    ## Compute covariance components
    

    

    
    if(! is.null(Z)){
      qr_decomp <- correct_int$qr
      R <- qr.R(qr_decomp)
      R_inv <- backsolve(R, diag(ncol(R)))
      const1 <-  R_inv %*% t(R_inv)
      C22 = N* (const1 %*% final_C[(1:(dim_Z + 1)), 1:(dim_Z + 1)] %*% const1)[1, 1]
      
      qr_decomp <- incorrect_int$qr
      R <- qr.R(qr_decomp)
      R_inv <- backsolve(R, diag(ncol(R)))
      const2 <- R_inv %*% t(R_inv)
      C33 = final_C[(dim_Z + 2):dim(final_C)[2],  (dim_Z + 2):dim(final_C)[2]]
      C33 = (const2 %*% C33 %*% const2 )[2, 2] * N
      
      
      C23 = final_C[1:(dim_Z + 1),  (dim_Z + 2):dim(final_C)[2]]
      C23 = ( (const1) %*% C23 %*% const2 )[1, 2] * N
      
    } else{
      const1 = solve(t(XZ_tilde) %*% XZ_tilde)
      C22 = N*  final_C[1, 1] * const1^2
      
      C33 = final_C[(dim_Z + 2):dim(final_C)[2],  (dim_Z + 2):dim(final_C)[2]]
      design = t(cbind(rep(1, N), XZ)) %*% cbind(rep(1, N), XZ)
      C33 = (solve(design) %*% C33 %*% solve(design) )[2, 2] * N
      
      
      C23 = final_C[1:(dim_Z + 1),  (dim_Z + 2):dim(final_C)[2]]
      C23 = ((const1) %*% matrix(C23, 1, 2) %*% solve(design) )[1, 2] * N
      
      
    }
    


    # Compute C12, C13
    
    C12 = overlap_ratio * C23 * N / N_ext
    C13 = overlap_ratio * C33 * N / N_ext
    
    # Covariance matrix of (alpha_ext - alpha_int, beta_int - beta)
    result_cov = rbind(c(C11 - 2 * C13+ C33, C12 - C23), c(C12 - C23, C22))
    # Compute the final calibrated estimator
    beta_cal = beta_int +  (C23 - C12) / ( C11 + C33 - 2 * C13) * (alpha_ext - alpha_int)
    beta_cal_var = (C22 - (C23 - C12)^2 / (C11 + C33 - 2 * C13) ) / N

    # Compute internal false model by sampling one sibling from each family
    
    sub_ind = sample_indices(F_ind, max(F_ind))
    
    incorrect_int = lm( Y[sub_ind] ~  XZ[sub_ind, ])
    alpha_int_single = summary(incorrect_int)$coefficients[2,1]
    alpha_int_var_single = summary(incorrect_int)$coefficients[2,2]^2
    
    return(list(beta_cal  = beta_cal, beta_cal_var = (beta_cal_var), beta_int = beta_int,
    beta_int_var = (C22/N), alpha_int = alpha_int, alpha_int_var = (C33/N),
    alpha_int_single = alpha_int_single, alpha_int_var_single = alpha_int_var_single) )
  }
    
    
    )
}
            

format_data2 <- function(pheno, F_ind, Z = NULL){

  F_ind = match(F_ind, unique(F_ind))
  
  order_indices <- order(F_ind)
  pheno <- pheno[order_indices]
  Z <- Z[order_indices, ]
  F_ind <- F_ind[order_indices]
  
  
  size_dic <- as.integer(table(F_ind))
  F_size <- size_dic[F_ind]
  
  # Reshape pheno and Z into a two-row format by family index
  pheno_matrix <- matrix(pheno, nrow = 2)
  Z_matrix <- array(Z, dim = c(2, nrow(Z) / 2, ncol(Z)))
  
  # Compute the differences along the two-row pairs
  pheno_diff <- pheno_matrix[1, ] - pheno_matrix[2, ]
  Z_diff <- Z_matrix[1, , ] - Z_matrix[2, , ]
  
  # Return results as a list

  dim_Z = dim(Z)
  
  return(list(Y = pheno, Y1 = pheno_matrix[1, ], Y2 = pheno_matrix[2, ], Y_diff = pheno_diff, F_ind = F_ind, Z = Z, 
              Z1 = Z_matrix[1, , ], Z2 = Z_matrix[2, , ], Z_diff = Z_diff, order = order_indices,
              size_dic = size_dic, F_size = F_size, K = max(F_ind),  N = length(F_ind), dim_Z = dim_Z))
}





            

calibrated_estimator2 <- function(X, data, alpha_ext, alpha_ext_var, N_ext,
                                 overlap_ratio = 0){
  with(data, {
    # Adjust external variance
    C11 = (alpha_ext_var / N_ext) * N
    
    # Compute correct model
    
    X = X[order]
    X = matrix(X, N, 1)

    X_matrix <- matrix(X, nrow = 2)
    X_diff = X_matrix[1, ] - X_matrix[2, ]
    
    X_diff = matrix(X_diff, N/2, 1)
    X1 = matrix(X_matrix[1, ], N/2, 1)
    
    XZ = cbind(X, Z)
    XZ_diff = cbind(X_diff, Z_diff)    
    
    correct_int = lm(Y_diff ~ -1 +  XZ_diff)
    beta_int = summary(correct_int)$coefficients[1, 1]
    beta_int_var = summary(correct_int)$coefficients[1, 2]^2

    
    # Compute the incorrect model for internal data
    
    incorrect_int = lm( Y ~  XZ)
    alpha_int = summary(incorrect_int)$coefficients[2, 1]
    


    # Compute internal false model by sampling one sibling from each family
    
    
    incorrect_int = lm( Y1 ~ cbind(X1, Z1) )
    alpha_int_single = summary(incorrect_int)$coefficients[2,1]
    alpha_int_var_single = summary(incorrect_int)$coefficients[2,2]^2
    
    return(list( beta_int = beta_int, beta_int_var = beta_int_var, alpha_int = alpha_int,
    alpha_int_single = alpha_int_single, alpha_int_var_single = alpha_int_var_single, N = N) )
  }
    
    
    )
}


