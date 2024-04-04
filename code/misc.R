# @param Y: a vector contains response variable values
# @param X1: a vector contains transmitted allele values
compute_sumstat <- function(Y, X1, family = "gaussian"){
  n <- length(Y)
  if (family == "gaussian"){
    fit <- lm(Y ~ X1)
  }

  if (family == "binomial"){
    fit <- glm(Y ~ X1, family = "binomial")
  }
  fit_summary = summary(fit)
  bhat = fit_summary$coefficients[2, "Estimate"]
  std = fit_summary$coefficients[2, "Std. Error"]
  resid <- fit$residuals
  return(list(bhat = bhat, std = std, resid = resid))
}
