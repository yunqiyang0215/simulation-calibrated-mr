

### Helper functions to simulate trio data.
# @par.allele: allele value of a parent
sample_allele <- function(par.allele){
  if (par.allele == 0) val = 0
  if (par.allele == 2) val = 1
  if (par.allele == 1) val = sample(c(0, 1), 1)
  return(val)
}

# @par1: genotyp of parent 1
# @par2: genotype of parent 2
# @return: transmitted and non-transmitted alleles of offspring
sim_offspring <- function(par1, par2){
  p = length(par1)
  offspring.t = rep(NA, p)
  offspring.nt = rep(NA, p)
  for (i in 1:p){
    a1 <- sample_allele(par1[i])
    a2 <- sample_allele(par2[i])
    t <- a1 + a2
    nt <- par1[i] + par2[i] - t
    offspring.t[i] <- t
    offspring.nt[i] <- nt
  }
  return(list(offspring.t = offspring.t, offspring.nt = offspring.nt))
}
