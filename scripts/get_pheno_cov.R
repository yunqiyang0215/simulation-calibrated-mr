# Description: get phenotype and covariates

library(data.table)
library(dplyr)

# Get traits for MR
setwd("/home/yunqiyang/calibrated_mr/real_data_analysis")
input.file <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                        "9-jan-2024","ukb676826.csv.gz")

dat <- fread(input.file,sep = ",",header = TRUE,verbose = FALSE,
             showProgress = FALSE,colClasses = "character")
class(dat) <- "data.frame"

cols <- c("eid", "4079-0.0", "4080-0.0", "21001-0.0", "2443-0.0", "6138-0.0", "1239-0.0")
col_names <- c("id", "dbp", "sbp", "bmi", "diabetes", "education_attainment", "smoking")
dat <- dat[, cols]
colnames(dat) <- col_names


# Get covariates
input.file <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                        "12-may-2023","ukb676949.csv.gz")

dat1 <- fread(input.file,sep = ",",header = TRUE,verbose = FALSE,
             showProgress = FALSE,colClasses = "character")
class(dat1) <- "data.frame"

cols <- c("eid","31-0.0","54-0.0","21022-0.0","22006-0.0",
          "22001-0.0","22000-0.0","22005-0.0",paste0("22009-0.",1:20),
          "22021-0.0", "22027-0.0")

col_names <- c("id","sex","assessment_centre","age","ethnic_genetic",
               "sex_genetic","genotype_measurement_batch","missingness",
               paste0("pc_genetic",1:20), "kinship_genetic", "outliers")

covar     <- dat1[,cols]
names(covar) <- col_names


# Merge two data
common_ids <- intersect(dat$id, covar$id)
dat_common <- dat[dat$id %in% common_ids, ]
covar_common <- covar[covar$id %in% common_ids, ]

# Sort data by IDs
dat_common <- dat_common[order(dat_common$id), ]
covar_common <- covar_common[order(covar_common$id), ]
sum(dat_common$id != covar_common$id)
dat = cbind(dat_common, covar_common[,c(2:ncol(covar_common))])

# Remove rows with mismatches between self-reported and genetic sex
dat <- dat %>% filter(sex == sex_genetic)
cat(sprintf("After removing sex mismatches, %d rows remain.\n",nrow(dat)))

# Remove any samples that are not marked as being "White British".
# This step should filter out 92887 rows.
dat = dat %>% filter(ethnic_genetic == 1)
cat(sprintf("After removing non White British, %d rows remain.\n",nrow(dat)))

# Remove "missingness" and "heterozygosity" outliers as defined by UK
# Biobank. Note that this step will remove any samples in which 
# the "missingness" column is greater than 5%.
dat <- dat %>% filter(outliers != 1)
cat(sprintf("After removing outliers, %d rows remain.\n",nrow(dat)))

saveRDS(dat, "./data/pheno_cov.rds")


