
# @param Y: a vector contains response variable values
# @param X1: a vector contains transmitted allele values
compute_sumstat <- function(Y, X1){
  n <- length(Y)
  mod <- lm(Y ~ X1)
  resid <- mod$residuals
  bhat <- mod$coef[2]
  return(list(bhat = bhat, resid = resid))
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






