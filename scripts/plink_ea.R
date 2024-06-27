library(dplyr)
dat = readRDS("/home/yunqiyang/calibrated_mr/real_data_analysis/data/pheno_cov.rds")

# Get GWAS samples. Remove relatedness
dat <- dat %>% filter(kinship_genetic == 0)
cat(sprintf(paste("After removing relatedness individuals based on kinship,",
                  "%d rows remain.\n"),nrow(dat)))

out.id.file <- "/scratch/yunqiyang/ukb_gwas_sumstat/education_yrs/id.txt"
out.pheno.file <- "/scratch/yunqiyang/ukb_gwas_sumstat/education_yrs/pheno.txt"
out.covar.file <- "/scratch/yunqiyang/ukb_gwas_sumstat/education_yrs/covar.txt"

# process EA 
dat = dat[!(dat$education_attainment %in% c(6, -7, -3)), ] 
dat$education_attainment = as.numeric(dat$education_attainment)
dat = dat[!is.na(dat$education_attainment), ] 

# map EA to years
dat$education_attainment[which(dat$education_attainment == 1)] = 17
dat$education_attainment[which(dat$education_attainment == 2)] = 14
dat$education_attainment[which(dat$education_attainment == 3)] = 11
dat$education_attainment[which(dat$education_attainment == 4)] = 11
dat$education_attainment[which(dat$education_attainment == 5)] = 15

IID = dat$id
FID = dat$id

# prepare id file
id.table = cbind(FID, IID)
write.table(id.table, file= out.id.file, quote=FALSE, row.names = FALSE, col.names = FALSE)

# prepare pheno file
pheno.table = data.frame(cbind(FID, IID, dat$education_attainment))
colnames(pheno.table) = c("FID", "IID", "education_yrs")
pheno.table$education_yrs = scale(as.numeric(pheno.table$education_yrs))
write.table(pheno.table, file= out.pheno.file, quote=FALSE, row.names = FALSE)

# prepare covariates
cols = which(colnames(dat) %in% c("age",paste0("pc_genetic", 1:10)))
Z = scale(apply(dat[, cols], 2, as.numeric))
Z = as.data.frame(Z)
Z$age2 = Z$age^2
Z$sex = factor(dat$sex)
cov.table = cbind(FID, IID, Z)
write.table(cov.table, file= out.covar.file, quote=FALSE, row.names = FALSE)


