# @param Y: a vector contains response variable values
# @param X1: a vector contains transmitted allele values
compute_sumstat <- function(Y, X1){
  n <- length(Y)
  fit <- lm(Y ~ X1)
  fit_summary = summary(fit)
  
  bhat = fit_summary$coefficients[2, "Estimate"]
  std = fit_summary$coefficients[2, "Std. Error"]
  resid <- fit$residuals
  return(list(bhat = bhat, std = std, resid = resid))
}


# @param Y: a vector contains response variable values
# @param X1: a vector contains transmitted allele values
# @param X2: a vector contains non-transmitted allele values
compute_adj_sumstat <- function(Y, X1, X2){
  n <- length(Y)
  mod <- lm(Y ~ X1 + X2)
  resid <- mod$residuals
  bhat <- c(mod$coef[2], mod$coef[3])
  return(list(bhat = bhat, resid = resid))
}


# This is to assume internal data is a subset of external full data.
# And we have access to internal data.
# @param dat_int: n by 3 matrix of internal full data, (Y1, X1, X2)
# @param N: sample size of external data
# @param resid: residuals from omitted variable model
# @param resid_adj: residuals from correct model.
# @param coef: coefficient from unadjusted model on internal data
# @param coef_adj: coefficient from adjusted model, (\hat b1, \hat b2)
# @param coef_ext: coefficient from unadjusted model on external data
calibrated_estimator <- function(dat_int, N, resid, resid_adj, coef, coef_adj, coef_ext){
  n <- nrow(dat_int)
  rho <- n/N

  X <- cbind(rep(1,n), dat_int[,2], dat_int[,3])
  X_tilde <- cbind(rep(1,n), dat_int[,2])
  sigma_hat_sq <- sum(resid_adj^2)/(n-3)
  sigma_star_hat_sq <- sum(resid^2)/(n-2)

  # this returns (X_i*X_i^T) which is a 3*3 matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # this returns (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix as the i-th column of out.prod2
  out.prod2 <- apply(X_tilde,1,function(x){outer(x,x)})
  D_1_hat <- - (1/(n*sigma_hat_sq))*matrix(rowSums(out.prod1),3,3)
  D_2_hat <- - (1/(n*sigma_star_hat_sq))*matrix(rowSums(out.prod2),2,2)
  C_11_hat <- (1/n)*(1/sigma_hat_sq^2)*matrix(rowSums(out.prod1%*%diag(resid_adj)^2),3,3)
  C_22_hat <- (1/n)*(1/sigma_star_hat_sq^2)*matrix(rowSums(out.prod2%*%diag(resid)^2),2,2)
  out.prod3 <- out.prod1[-c(7,8,9),]
  C_12_hat <- (1/n)*(1/(sigma_hat_sq*sigma_star_hat_sq))*
    matrix(rowSums(out.prod3%*%diag(resid*resid_adj)),3,2,byrow = "F")

  a <- solve(D_1_hat) %*% C_11_hat %*% solve(D_1_hat)
  b <- solve(D_1_hat) %*% C_12_hat %*% solve(D_2_hat)
  c <- solve(D_2_hat) %*% C_22_hat %*% solve(D_2_hat)
  V_1 <- cbind(a,b,rho*b)
  V_2 <- cbind(t(b),c,rho*c)
  V_3 <- cbind(rho*t(b),rho*c,rho*c)
  V <- rbind(V_1,V_2,V_3)
  A <- rbind(c(0,1,-1,0,0,0,0),c(0,0,0,0,1,0,-1))
  V_hat <- A%*%V%*%t(A)
  # calculate calibrated estimator
  tau_cal <- coef_adj[1] - coef_adj[2] - V_hat[1,2]/V_hat[2,2]*(coef - coef_ext)
  var_tau_cal <- (V_hat[1, 1] - V_hat[1,2]^2/V_hat[2,2])/n

  # calculate raw estimator
  tau_raw <- coef_adj[1] - coef_adj[2]
  var_tau_raw <- (V_hat[1, 1])/n
  return(list(tau_cal = tau_cal, tau_raw = tau_raw, var_tau_cal = var_tau_cal, var_tau_raw = var_tau_raw))
}


# @param data_int: 3 columns, Y, X1 and X2 where X1 = G and X2 = G_PA
# @param data_ext: 2 columns, Y, X1
# @param rho_shared = n_shared/N
# @param alpha_ext = summary statistic gwas
# @param sd_alpha_ext
calibrated_est_trio_2 <- function(data_int,rho_shared,alpha_ext,sd_alpha_ext){
  # Creating the internal dataset
  Y_int <- data_int[,1]
  X_1_int <- data_int[,2]
  X_2_int <- data_int[,3]
  #internal model
  mod_int <- lm(Y_int ~ X_1_int + X_2_int)
  #external model
  mod_unadjusted <- lm(Y_int ~ X_1_int)
  n <- nrow(data_int)
  # creating the X matrix for internal model
  X <- cbind(rep(1,n),X_1_int,X_2_int)
  # creating the X_tilde matrix for external data
  X_tilde <- cbind(rep(1,n),X_1_int)
  sigma_hat_sq <- sum(mod_int$residuals^2)/(n-3)
  sigma_star_hat_sq <- sum(mod_unadjusted$residuals^2)/(n-2)
  # this returns (X_i*X_i^T) which is a 3*3 matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # this returns (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix
  # as the i-th column of out.prod2
  #out.prod2 <- apply(X_tilde,1,function(x){outer(x,x)})
  out.prod2 <- out.prod1[c(1,2,4,5),]
  D_1_hat <- - (1/(n*sigma_hat_sq))*matrix(rowSums(out.prod1),3,3)
  D_2_hat <- - (1/(n*sigma_star_hat_sq))*matrix(rowSums(out.prod2),2,2)
  C_11_hat <- (1/n)*(1/sigma_hat_sq^2)*
    matrix(rowSums(out.prod1%*%diag((mod_int$residuals)^2)),3,3)
  C_22_hat <- (1/n)*(1/sigma_star_hat_sq^2)*
    matrix(rowSums(out.prod2%*%diag((mod_unadjusted$residuals)^2)),2,2)
  out.prod3 <- out.prod1[-c(7,8,9),]
  C_12_hat <- (1/n)*(1/(sigma_hat_sq*sigma_star_hat_sq))*
    matrix(rowSums(out.prod3%*%diag((mod_unadjusted$residuals)*(mod_int$residuals))),3,2,byrow = "F")
  D_1_inv <- solve(D_1_hat)
  D_2_inv <- solve(D_2_hat)
  a <- D_1_inv %*% C_11_hat %*% D_1_inv
  b <- D_1_inv %*% C_12_hat %*% D_2_inv
  c <- D_2_inv %*% C_22_hat %*% D_2_inv
  V_1 <- cbind(a,b,rho_shared*b)
  V_2 <- cbind(t(b),c,rho_shared*c)
  V_3 <- cbind(rho_shared*t(b),rho_shared*c,rho_shared*c)
  V <- rbind(V_1,V_2,V_3)
  V[7,7] <- n*sd_alpha_ext^2
  A <- rbind(c(0,1,0,0,0,0,0),c(0,0,0,0,1,0,-1))
  V_hat <- A%*%V%*%t(A)
  tau_raw <- as.vector(mod_int$coef)[2]
  tau_cal <-  tau_raw -
    V_hat[1,2]/V_hat[2,2]*(as.vector(mod_unadjusted$coef)[2] - alpha_ext)
  var_tau_cal <- V_hat[1,1] - V_hat[1,2]^2/V_hat[2,2]
  return (list("raw_est" = tau_raw, "raw_variance" = V_hat[1,1], "calibrated_est" = tau_cal,"cali_variance" = var_tau_cal))
}

# @param data_int: n by 3 matrix of internal full data, (Y1, X1, X2)
# @param N: sample size of external data
# @param alpha_ext: coefficient from unadjusted model on external data
# @param sd_alpha_ext: sd of alpha_ext from GWAS summary statistic
# @param rho_shared: proportion of shared samples between internal and external data
# @returns calibrated estimator for logistic regression in trio data
calibrated_estimator_logistic <- function(data_int, rho_shared, N, alpha_ext, sd_alpha_ext){
  Y_int <- data_int[,1]
  X_1_int <- data_int[,2]
  X_2_int <- data_int[,3]
  n <- nrow(data_int)
  ##internal model
  mod_int <- glm(Y_int ~ X_1_int + X_2_int, family = binomial)
  ##external model
  mod_unadjusted <- glm(Y_int ~ X_1_int, family = binomial)
  #coef_adjusted <- c(mod_int$coef[2], mod_int$coef[3])
  #coef <- mod_ext$coef[2]

  X <- cbind(rep(1,n), X_1_int, X_2_int)
  X_tilde <- cbind(rep(1,n), X_1_int)

  # this returns (X_i*X_i^T) which is a 3*3 matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # this returns (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix as the i-th column of out.prod2
  out.prod2 <- out.prod1[c(1,2,4,5),]
  D_1_hat <- - (1/n)*matrix(rowSums(out.prod1%*%diag(mod_int$fitted.values*(1 - mod_int$fitted.values))),3,3)
  D_2_hat <- - (1/n)*matrix(rowSums(out.prod2%*%diag(mod_unadjusted$fitted.values*(1 - mod_unadjusted$fitted.values))),2,2)
  C_11_hat <- (1/n)*matrix(rowSums(out.prod1%*%diag((data_int[,1] - mod_int$fitted.values)^2)),3,3)
  C_22_hat <- (1/n)*matrix(rowSums(out.prod2%*%diag((data_int[,1] - mod_unadjusted$fitted.values)^2)),2,2)
  out.prod3 <- out.prod1[-c(7,8,9),]
  C_12_hat <- (1/n)*matrix(rowSums(out.prod3%*%diag((data_int[,1] - mod_int$fitted.values)*
                                                      (data_int[,1] - mod_unadjusted$fitted.values))),3,2,byrow = "F")
  D_1_inv <- solve(D_1_hat)
  D_2_inv <- solve(D_2_hat)
  a <-  D_1_inv %*% C_11_hat %*% D_1_inv
  b <- D_1_inv %*% C_12_hat %*% D_2_inv
  c <- D_2_inv %*% C_22_hat %*% D_2_inv

  V_1 <- cbind(a,b,rep(0,3))
  V_2 <- cbind(t(b),c,rep(0,2))

  V <- rbind(V_1,V_2,rep(0,6))
  V[6,6] <- n*sd_alpha_ext^2
  A <- rbind(c(0,1,-1,0,0,0),c(0,0,0,0,1,-1))
  V_hat <- A%*%V%*%t(A)
  # calculate calibrated estimator
  tau_raw <- as.vector(mod_int$coef)[2] - as.vector(mod_int$coef)[3]
  tau_cal <-  tau_raw -
    V_hat[1,2]/V_hat[2,2]*(as.vector(mod_unadjusted$coef)[2] - alpha_ext)
  var_tau_cal <- (1/n)*(V_hat[1,1] - V_hat[1,2]^2/V_hat[2,2])
  return (list("raw_est" = tau_raw,"raw_variance" = V_hat[1,1]/n,"calibrated_est" = tau_cal,"cali_variance" = var_tau_cal))
}



