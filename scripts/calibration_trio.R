library(dplyr)
source("/home/yunqiyang/calibrated_mr/trio_data_cov.R")

args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]


if (length(args) > 1) {
  val_to_remove <- eval(parse(text = args[2]))
} else {
  val_to_remove <- NULL
}

if (pheno == "diabetes"){
  val_to_remove = c(-1, -3)
}

if (pheno == "education_yrs"){
  val_to_remove = c(6, -7, -3)
}


dat_int = readRDS("/home/yunqiyang/calibrated_mr/real_data_analysis/data/trio_geno_pheno.rds")
pheno.file = paste0("/scratch/yunqiyang/ukb_gwas_sumstat/gwas_", pheno, ".linear")
pheno.gwas = read.csv(pheno.file, sep = "\t")


var_names = colnames(dat_int)[grep("^rs", colnames(dat_int))]
ss = strsplit(var_names, split = "_")
snp_list = sapply(ss, function(x) x[1])
alleles = sapply(ss, function(x) x[2])

## processing dat_int
# rename SNP
colnames(dat_int) <- gsub("_G|_T|_A|_C", "", colnames(dat_int))
dat_int <- data.frame(lapply(dat_int, function(x) as.numeric(as.character(x))), stringsAsFactors = FALSE)
# remove NA
dat_int <- dat_int[!is.na(dat_int[, pheno]), ]

# val_to_remove
if (!is.null(val_to_remove)){
  dat_int = dat_int[!(dat_int[, pheno] %in% val_to_remove), ]
}

# create output
p <- length(snp_list)
sumstat = matrix(NA, nrow = p, ncol = 8)
colnames(sumstat) = c("variant", "allele.test", "cali", "cali.var", "raw", "raw.var", "int", "int.var")


dat_child = dat_int[dat_int$If.Child == 1, ]
dat_par = dat_int[dat_int$If.Child == 0, c("Individual.ID", "Family.ID", snp_list)]

# only keep families with two parents
dat_par <- dat_par %>%
  group_by(Family.ID) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  as.data.frame()

# get parents genotype sum
X.pa = c()
family_indx = unique(dat_par$Family.ID)
for (indx in family_indx){
  dat_sub = dat_par[dat_par$Family.ID == indx, snp_list]
  val = colSums(dat_sub, na.rm = TRUE)
  X.pa = rbind(X.pa, c(indx, val))
}
colnames(X.pa)[1] = "Family.ID"

dat = merge(dat_child, X.pa, by = "Family.ID", all = FALSE)

# Using the same covariates adjusted in GWAS. Y is scaled 
cols = which(colnames(dat) %in% c("age", paste0("pc_genetic", 1:10)))
Z = scale(apply(dat[, cols], 2, as.numeric))
Z = as.data.frame(Z)
Z$age2 = Z$age^2
Z$sex = dat$sex
Z = as.matrix(Z)

Y = scale(dat[, pheno])

N_ext = 2e5
for (i in 1:length(snp_list)){
  snp = snp_list[i]
  if (!(snp %in% pheno.gwas$ID)){
    next
  }
  indx = which(pheno.gwas$ID == snp)

  # check if the allele tested in external gwas the same as internal genotype dosage
  alpha_ext = pheno.gwas[indx, "BETA"]
  alpha_ext_sd = pheno.gwas[indx, "SE"]
 
  if (pheno.gwas[indx, "A1"] != alleles[i]){
    alpha_ext = -alpha_ext
  }
  
  snp_col = c(paste0(snp, ".x"), paste0(snp, ".y"))
  snp_columns <- which(colnames(dat) == snp_col[1] | colnames(dat) == snp_col[2])
  dat_exp = cbind(Y, dat[ , snp_columns])
  res <- calibrated_est_trio_2(dat_exp, Z = Z, rho_shared = 0, alpha_ext, alpha_ext_sd)
  
  # get internal fit 
  dat_complete = as.data.frame(cbind(Y, dat[ , snp_columns][,1], Z))
  colnames(dat_complete)[1] = "Y" 
  fit = lm(Y ~ ., data = dat_complete)
  int = coef(summary(fit))[2, 1]
  int.var = coef(summary(fit))[2,2]^2
  
  sumstat[i, ] = c(snp, alleles[i], res$calibrated_est, res$variance, res$`Raw Estimator`, res$raw_variance, int, int.var)
}

saveRDS(sumstat, paste0("/home/yunqiyang/calibrated_mr/real_data_analysis/result/trio_calibration_", pheno, ".rds"))



