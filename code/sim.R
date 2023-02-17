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
