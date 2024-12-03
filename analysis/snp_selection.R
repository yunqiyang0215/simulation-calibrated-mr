---
title: "MR SNP selection"
author: "Yunqi Yang"
date: "08/07/2023"
output: html_document
---

library(mr.raps)
library(GRAPPLE)
library(MendelianRandomization)

### 1. height on EA
sel.file <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/Height-giant14.csv"
exp.file <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/height-ukb.csv"
out.file <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/GWAS_EA_excl23andMe.txt"
plink_refdat <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/data_maf0.01_rs_ref/data_maf0.01_rs_ref"

pval = c(1e-2, 1e-5, 1e-8)
for (i in 1:length(pval)){

  data.list <- getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres = pval[i],
                        plink_exe = '/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/plink_mac_20230116/plink')
  write.table(data.list$data[, 1], paste0("/Users/nicholeyang/Downloads/simulation-calibrated-mr/data/snp_list/height_",
                                          pval[i], ".txt"), quote = FALSE, col.names  = FALSE, row.names = FALSE)
}




### 2. BMI on education attainment
sel.file <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/BMI-giant17eu.csv"
exp.file <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/BMI-ukb.csv"

exp.file <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/BMI-giant17eu.csv"
out.file <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/GWAS_EA_excl23andMe.txt"
plink_refdat <- "/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/data_maf0.01_rs_ref/data_maf0.01_rs_ref"

pval = c(1e-2, 1e-5, 1e-8)
for (i in 1:length(pval)){

  data.list <- getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres = pval[i],
                      plink_exe = '/Users/nicholeyang/Downloads/calibrated_estimator/simulation/real_data/plink_mac_20230116/plink')
  write.table(data.list$data[, 1], paste0("/Users/nicholeyang/Downloads/simulation-calibrated-mr/data/snp_list/bmi_",
                                          pval[i], ".txt"), quote = FALSE, col.names  = FALSE, row.names = FALSE)
}


