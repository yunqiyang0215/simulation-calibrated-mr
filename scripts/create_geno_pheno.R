# Description: create internal phenotype data
library(dplyr)

# create trio data 
# Step 1: filter the trio data and merge with phenotype data
trio = readRDS("/home/yunqiyang/calibrated_mr/ukb_family/trio_with_fam_structure.rds")
pheno = readRDS("/home/yunqiyang/calibrated_mr/real_data_analysis/data/pheno_cov.rds")

## Merge two data
trio = data.frame(trio)
common_ids <- intersect(trio$Individual.ID, pheno$id)
trio_common <- trio[trio$Individual.ID %in% common_ids, ]
pheno_common <- pheno[pheno$id %in% common_ids, ]

## Sort data by IDs
trio_common <- trio_common[order(trio_common$Individual.ID), ]
pheno_common <- pheno_common[order(pheno_common$id), ]
sum(trio_common$Individual.ID != pheno_common$id)
dat = cbind(trio_common, pheno_common[,c(2:ncol(pheno_common))])

## Keep families with exact 3 members
dat <- dat %>%
  group_by(Family.ID) %>%
  filter(n() == 3) %>%
  ungroup()


# step 2: merge with genotype data
geno = readRDS("/scratch/yunqiyang/geno_ukb/bmi/trio_geno.rds")
dat2 = merge(dat, geno, by.x = "Individual.ID", by.y = "IID")

# create education yrs 
dat2$education_yrs = dat2$education_attainment
dat2$education_yrs[which(dat2$education_attainment == 1)] = 17
dat2$education_yrs[which(dat2$education_attainment == 2)] = 14
dat2$education_yrs[which(dat2$education_attainment == 3)] = 11
dat2$education_yrs[which(dat2$education_attainment == 4)] = 11
dat2$education_yrs[which(dat2$education_attainment == 5)] = 15
saveRDS(dat2, "/home/yunqiyang/calibrated_mr/real_data_analysis/data/trio_geno_pheno.rds")





# Step 1: filter the sibling data and merge with phenotype data
sib = read.csv("/home/yunqiyang/calibrated_mr/ukb_family/fam_sib.csv")
pheno = readRDS("/home/yunqiyang/calibrated_mr/real_data_analysis/data/pheno_cov.rds")

## Merge two data
common_ids <- intersect(sib$Individual, pheno$id)
sib_common <- sib[sib$Individual %in% common_ids, ]
pheno_common <- pheno[pheno$id %in% common_ids, ]

## Sort data by IDs
sib_common <- sib_common[order(sib_common$Individual), ]
pheno_common <- pheno_common[order(pheno_common$id), ]
sum(sib_common$id != pheno_common$id)
dat = cbind(sib_common, pheno_common[,c(2:ncol(pheno_common))])

## Keep families with exact two siblings
dat = dat[dat$Size == 2, ]
dat <- dat %>%
  group_by(Family.ID) %>%
  filter(n() > 1) %>%
  ungroup()

# step 2: merge with genotype data
geno = readRDS("/scratch/yunqiyang/geno_ukb/bmi/sib_geno.rds")
dat2 = merge(dat, geno, by.x = "Individual", by.y = "IID")
# create education yrs
dat2$education_yrs = dat2$education_attainment
dat2$education_yrs[which(dat2$education_attainment == 1)] = 17
dat2$education_yrs[which(dat2$education_attainment == 2)] = 14
dat2$education_yrs[which(dat2$education_attainment == 3)] = 11
dat2$education_yrs[which(dat2$education_attainment == 4)] = 11
dat2$education_yrs[which(dat2$education_attainment == 5)] = 15
saveRDS(dat2, "/home/yunqiyang/calibrated_mr/real_data_analysis/data/sib_geno_pheno.rds")




