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

##' Compute the calibrated estimator
##' @param Y: a vector that contains response variable values
##' @param X: a vector that contains transmitted allele values
##' @param Z: a matrix that contains other covariates
##' @param F_ind: family index (from 1 to K)
##' @param alpha_ext: estimate of alpha from external data
##' @param alpha_ext_var: variance of alpha from external data
##' @param N_ext: number of samples in external data
##' @param overlap_ratio: proportion of internal data from external data
##' @param Z: other covariates (the number of rows should be internal data size)
##' 
calibrated_estimator <- function(Y, X,F_ind, alpha_ext, alpha_ext_var, N_ext,  
                                 overlap_ratio = 0, Z = NULL){
  
  # Correct model Y = beta_0 + beta_1 X + f + epsilon
  # Incorrect model: Y = alpha_0 + alpha_1 X + epsilo
  
  # Number of families
  K = max(F_ind)
  # Number of individuals
  N = length(X)
  
  size_dic = rep(0, K) # k-th entry represents number of individuals in family k
  F_size = rep(0, N) # n-th entry represents the number of people in n-th person's family
  for(i in 1:N){
    size_dic[F_ind[i]] = size_dic[F_ind[i]]  + 1
  }
  
  for(i in 1:N){
    F_size[i] = size_dic[F_ind[i]]
  }
  
  
  # Adjust external variance
  C11 = (alpha_ext_var / N_ext) * N
  
  # Compute the correct model for internal data
  
  # Y_tilde = Y - Y_bar
  # X_tilde = X - X_bar
  # Y_tilde = \beta_1 X_tilde + epsilon
  
  if(! is.null(Z)){
    dim_Z = length(Z) / N
  } else{
    dim_Z = 0
  }
  Y_tilde = Y
  
  X = matrix(X, N, 1)
  XZ = cbind(X, Z)
  XZ_tilde = XZ
  
  if(dim_Z != 0){
    for(i in 1:K){
      Y_tilde[ F_ind == i ] = Y_tilde[ F_ind == i ]  - mean(Y_tilde[ F_ind == i ] )
      XZ_tilde[ F_ind == i, ] = sweep(XZ_tilde[ F_ind == i, ], 2, 
                                      colMeans(XZ_tilde[ F_ind == i, ]), FUN="-")     }
  } else{
    for(i in 1:K){
      Y_tilde[ F_ind == i ] = Y_tilde[ F_ind == i ]  - mean(Y_tilde[ F_ind == i ] )
      XZ_tilde[ F_ind == i] = XZ_tilde[ F_ind == i]  - mean(XZ_tilde[ F_ind == i] )
    }
  }
  
  XZ = matrix(XZ, ncol = 1 + dim_Z)
  XZ_tilde = matrix(XZ_tilde, ncol = 1 + dim_Z)
  
  
  correct_int = lm(Y_tilde ~ -1 +  XZ_tilde)
  beta_int = summary(correct_int)$coefficients[1]
  
  # Compute the incorrect model for internal data
  
  incorrect_int = lm( Y ~  XZ)
  alpha_int = summary(incorrect_int)$coefficients[2]
  # Compute the C_{22}, C_{33}, C_{23}
  
  M = max(F_size)
  
  Cs = list()
  
  # Consider each family size separately
  for(m in 2:M){
    
    # Only look at family of size m
    XZ_sub = XZ[F_size == m, ]
    resid_sub = resid(summary(incorrect_int))[F_size == m]
    
    XZ_tilde_sub = XZ_tilde[F_size == m, ]
    resid_tilde_sub = resid(summary(correct_int))[F_size == m]
    
    
    ind_sub = F_ind[F_size == m]
    
    C = cbind(XZ_tilde_sub * resid_tilde_sub, resid_sub, XZ_sub * resid_sub)
    
    C = t(sapply(unique(ind_sub), 
                 function(x) colSums(C[which(ind_sub == x), ])))
    
    
    Cs[[m]] = cov(C)
    
  }  
  
  final_C = matrix(0, 2 * (dim(XZ)[2]) + 1, 2 * (dim(XZ)[2]) + 1)
  
  for(i in 1:K){
    final_C = final_C + Cs[[size_dic[i]]]
  }
  
  
  ## Compute covariance components
  
  if(! is.null(Z)){
    const1 = solve(t(XZ_tilde) %*% XZ_tilde)
    C22 = N* (const1 %*% final_C[(1:(dim_Z + 1)), 1:(dim_Z + 1)] %*% const1)[1, 1]
    
    C33 = final_C[(dim_Z + 2):dim(final_C)[2],  (dim_Z + 2):dim(final_C)[2]]
    design = t(cbind(rep(1, N), XZ)) %*% cbind(rep(1, N), XZ)
    C33 = (solve(design) %*% C33 %*% solve(design) )[2, 2] * N
    
    
    C23 = final_C[1:(dim_Z + 1),  (dim_Z + 2):dim(final_C)[2]]
    C23 = ( (const1) %*% C23 %*% solve(design) )[1, 2] * N
    
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
  return(list(beta_cal  = beta_cal, beta_cal_var = (beta_cal_var), beta_int = beta_int, 
              beta_int_var = (C22/N)))
}







# Example
# n_trial = 1000
# K_ext = 50000
# K_int = 200
#
# beta_int = rep(0, n_trial)
# beta_cal = rep(0, n_trial)
#
# beta_cal_sd = rep(0, n_trial)
# sigma_f = 1
# for(i in 1:n_trial){
#   external_data = generate_data(K = K_ext, sigma_f = sigma_f, seed = i, eps_cor = 0)
#   model = summary(lm(external_data$Y[1:K_ext] ~  external_data$X[1:K_ext]))
#
#   ext_estimate = model$coefficients[2, 1]
#   ext_var = model$coefficients[2, 2]^2
#   internal_data =  generate_data(K = K_int, sigma_f = sigma_f, seed = n_trial+ i, eps_cor = 0)
#   result = calibrated_estimator(internal_data$Y, internal_data$X, internal_data$F_ind, ext_estimate
#                                 , ext_var)
#   beta_int[i] = result$beta_int
#   beta_cal[i] = result$beta_cal
#   beta_cal_sd[i] = result$beta_cal_sd
# }
#
# # MSE
#
# sum( (beta_int - 1)^2)
# sum( (beta_cal - 1)^2)
#
# # Variance reduction
# var(beta_int)
# var(beta_cal)
#
# 1 - var(beta_cal)/var(beta_int)
#
# # Computed sd vs true sd
#
# sd(beta_cal)
# mean(beta_cal_sd)
#
# # Density plot
# plot(density(beta_int), ylim = c(0,20))
# lines(density(beta_cal), col = "red")

