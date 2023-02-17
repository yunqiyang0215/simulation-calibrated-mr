
# Simulation 1: check the derivation of variance reduction in Gaussian case (x1, x2 joint Gaussian)
rs <- c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)

source("/project2/mstephens/yunqiyang/calibrated_mr/sim_gaussian/sim_func.R")
library(mvtnorm)

t1 <- proc.time()

seeds <- c(1:1e5)
n = 5e3
N = 5e4
b = c(1,3,0.1)
sigma = 1
res.all <- list()


for (j in 1:length(rs)){
  r <- rs[j]
  res <- matrix(NA, ncol = 2, nrow = length(seeds))
  colnames(res) <- c("gamma.tilde", "gamma.raw")

  for (i in 1:length(seeds)){
    set.seed(i)
    # simulate data
    datN <- sim_dat(N, b, r = r, sigma = sigma)
    datn <- sim_dat(n, b, r = r, sigma = sigma)

    # fit 3 models
    fit.true <- lm(y ~ x1 + x2, dat = data.frame(datn))
    fitN <- lm(y ~ x1, dat = data.frame(datN))
    fitn <- lm(y ~ x1, dat = data.frame(datn))

    bhat <- coef(fit.true)[-1]
    ahat.N <- coef(fitN)[-1]
    ahat.n <- coef(fitn)[-1]

    res[i, 1] <- calibrated_estimator(bhat, ahat.N, ahat.n, rho = n/N, b[3], r = r, sigma = sigma)
    res[i, 2] <- bhat[1] - bhat[2] # \hat b1 - \hat b2
  }
  res.all[[j]] <- res
}

t2 <- proc.time()
print(t2 - t1)

saveRDS(res.all, "/project2/mstephens/yunqiyang/calibrated_mr/sim_gaussian/gaussian_r.rds")

