# Function to simulate X, assuming X is standardized, var(x1) = var(x2) = 1.
# If dist == "uniform" or "binomial", then X is first simulated from Gaussian, and then
# transformed to uniform/binomial distribution.
# @param n: sample size
# @param r: correlation between x1 and x2
# @param dist: the distribution of x1 and x2
sim_X <- function(n, r, dist = "gaussian"){
  cov = matrix(c(1,r,r,1), nrow = 2, ncol = 2)
  X <- rmvnorm(n, mean = rep(0, nrow(cov)), sigma = cov)
  if (dist == "uniform"){
    X <- pnorm(X)
  }
  if (dist == "binomial"){
    X[X > 0] = 1
    X[X <= 0] = 0
  }
  return(X)
}


# Function to simulte y and X from the true model,
# y = b0 + b1*x1 + b2*x2 + e, e ~ N(0, sigma^2).
# @param n: sample size
# @param b: vector of length 3, b = (b0, b1, b2)
# @param r: correlation between x1 and x2
# @param dist: the distribution of x1 and x2
sim_dat <- function(n, b, r, sigma, dist = "gaussian"){
  # Simulate X
  X <- sim_X(n, r, dist)
  e <- rnorm(n, sd = sigma)
  y <- cbind(rep(1, n), X) %*% b + e # add intercept
  dat <- cbind(y, X)
  colnames(dat) <- c("y", "x1", "x2")
  return(dat)
}


# Function to calculate the calibrated estimator using equation from
# section 3.2 of the write-up. v12 and v22 are the true values, instead of plug-in estimates.
# @param bhat: a vector of length 2, containing \hat b1 and \hat b2.
# @param ahat.N: the estimate from misspecified model on larger dataset.
# @param ahat.n: the estimate from misspecified model on small dataset.
# @param rho: the ratio of n/N
# @param b2: true effect size of x2
# @param sigma: sd of random error
# @param r: correlation between x1 and x2. r can be either the true r if
# simulated from Gaussian, or estimated r if x1 and x2 are transformed to uniform
calibrated_estimator <- function(bhat, ahat.N, ahat.n, rho, b2, r, sigma){
  v12 = sigma^2*(1-rho)
  v22 = sigma^2 + b2^2*(1-r^2)
  gamma.tilde <- bhat[1] - bhat[2] - v12*(1/v22)*(ahat.n - ahat.N)
  return(gamma.tilde)
}

# @param data_valid: 3 columns, Y, X1 and X2
# @param data_ext: 2 columns, Y, X1
calibrated_est <- function(data_valid,data_ext){
  # Creating the internal dataset
  Y_valid <- data_valid[,1]
  X_1_valid <- data_valid[,2]
  X_2_valid <- data_valid[,3]
  # Creating the external dataset
  Y_ext <- data_ext[,1]
  X_1_ext <- data_ext[,2]
  #internal model
  mod_int <- lm(Y_valid ~ X_1_valid + X_2_valid)
  #external model
  mod_ext <- lm(Y_ext[c(1:n)] ~ X_1_ext[c(1:n)])
  mod_ext1 <- lm(Y_ext ~ X_1_ext)
  n <- nrow(data_valid)
  N <- nrow(data_ext)
  rho <- n/N
  # creating the X matrix for internal model
  X <- cbind(rep(1,n),X_1_valid,X_2_valid)
  # creating the X_tilde matrix for external data
  X_tilde <- cbind(rep(1,n),X_1_ext[c(1:n)])
  sigma_hat_sq <- sum(mod_int$residuals^2)/(n-3)
  sigma_star_hat_sq <- sum(mod_ext$residuals^2)/(n-2)
  # this returns (X_i*X_i^T) which is a 3*3 matrix as the i-th column of out.prod1
  out.prod1 <- apply(X,1,function(x){outer(x,x)})
  # this returns (\tilde{X}_i*\tilde{X}_i^T) which is a 2*2 matrix as the i-th column of out.prod2
  out.prod2 <- apply(X_tilde,1,function(x){outer(x,x)})
  D_1_hat <- - (1/(n*sigma_hat_sq))*matrix(rowSums(out.prod1),3,3)
  D_2_hat <- - (1/(n*sigma_star_hat_sq))*matrix(rowSums(out.prod2),2,2)
  C_11_hat <- (1/n)*(1/sigma_hat_sq^2)*matrix(rowSums(out.prod1%*%diag((mod_int$residuals)^2)),3,3)
  C_22_hat <- (1/n)*(1/sigma_star_hat_sq^2)*matrix(rowSums(out.prod2%*%diag((mod_ext$residuals)^2)),2,2)
  out.prod3 <- out.prod1[-c(7,8,9),]
  C_12_hat <- (1/n)*(1/(sigma_hat_sq*sigma_star_hat_sq))*
    matrix(rowSums(out.prod3%*%diag((mod_ext$residuals)*(mod_int$residuals))),3,2,byrow = "F")
  a <- solve(D_1_hat) %*% C_11_hat %*% solve(D_1_hat)
  b <- solve(D_1_hat) %*% C_12_hat %*% solve(D_2_hat)
  c <- solve(D_2_hat) %*% C_22_hat %*% solve(D_2_hat)
  V_1 <- cbind(a,b,rho*b)
  V_2 <- cbind(t(b),c,rho*c)
  V_3 <- cbind(rho*t(b),rho*c,rho*c)
  V <- rbind(V_1,V_2,V_3)
  A <- rbind(c(0,1,-1,0,0,0,0),c(0,0,0,0,1,0,-1))
  V_hat <- A%*%V%*%t(A)
  tau_cal <- as.vector(mod_int$coef)[2] - as.vector(mod_int$coef)[3] -   V_hat[1,2]/V_hat[2,2]*(as.vector(mod_ext$coef)[2] - as.vector(mod_ext1$coef)[2])
  # calculate raw estimator
  bhat <- coef(mod_int)[-1]
  tau_raw <- bhat[1] - bhat[2]
  var_tau_cal <- (V_hat[1, 1] - V_hat[1,2]^2/V_hat[2,2])/n
  var_tau_raw <- (V_hat[1, 1])/n
  return(list(tau_cal = tau_cal, tau_raw = tau_raw, var_tau_cal = var_tau_cal, var_tau_raw = var_tau_raw))
}

