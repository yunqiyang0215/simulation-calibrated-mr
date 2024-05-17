library(dplyr)
source("/home/yunqiyang/calibrated_mr/Sibling_Calibrated_Estimator.R")
source("/home/yunqiyang/calibrated_mr/Sibling_Logistic_Estimator.R")


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



dat_int = readRDS("/home/yunqiyang/calibrated_mr/real_data_analysis/data/sib_geno_pheno.rds")
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

# keep families with 2 siblings
dat_int <- dat_int %>%
  group_by(Family.ID) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  as.data.frame()


# Using the same covariates adjusted in GWAS. Y is scaled 
cols = which(colnames(dat_int) %in% c("age", paste0("pc_genetic", 1:10)))
Z = scale(apply(dat_int[, cols], 2, as.numeric))
Z = as.data.frame(Z)
Z$age2 = Z$age^2
Z$sex = dat_int$sex
Y = scale(dat_int[, pheno])

# Conduct calibration
p <- length(snp_list)
sumstat = matrix(NA, nrow = p, ncol = 8)
colnames(sumstat) = c("variant", "allele.test", "cali", "cali.var", "raw", "raw.var", "int", "int.var")


N_ext = 2e5

for (i in 1:length(snp_list)){
  snp = snp_list[i]
  if (!(snp %in% pheno.gwas$ID)){
    next
  }
  indx = which(pheno.gwas$ID == snp)
  
  # check if the allele tested in external gwas the same as internal genotype dosage
  alpha_ext = pheno.gwas[indx, "BETA"]
  if (pheno.gwas[indx, "A1"] != alleles[i]){
    alpha_ext = -alpha_ext
  }
  alpha_ext_var = pheno.gwas[indx, "SE"]^2
  
  F_ind = dat_int$Family.ID
  X = dat_int[, snp]
  res = calibrated_estimator(Y, X, F_ind, alpha_ext, alpha_ext_var, N_ext, Z = Z)
  sumstat[i, ] = c(snp, alleles[i], res$beta_cal, res$beta_cal_var, res$beta_int, 
                   res$beta_int_var, res$alpha_int, res$alpha_int_var)
}

saveRDS(sumstat, paste0("/home/yunqiyang/calibrated_mr/real_data_analysis/result/sib_calibration_", pheno, ".rds"))



