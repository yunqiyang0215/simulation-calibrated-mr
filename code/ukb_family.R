library(data.table)
dat <- fread("/gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb_rel_a27386_s488130.dat")

# Relationship inference criteria from King paper:
# https://www.chen.kingrelatedness.com/publications/pdf/BI26_2867.pdf

# Monozygotic twin: kinship > 1/2^(3/2), IBD0 < 0.1
# Trio & full sib kinship range: 1/2^(5/2) ~ 1/2^(3/2)
# Trio IBD0: < 0.1
# Full sib IBD0: 0.1 ~ 0.365
# Unrelated: kinship close to 0 and IBD0: > 1- 1/2^(5/2), close to 1.
# UKB paper use IBD0 < 0.0012 as parent-offspring pair


rk1 = 1/(2^(5/2))
rk2 = 1/(2^(3/2))


r1 = 0.0012
r2 = 0.365

# check unique pairs: 107161. Data total row = 107161
dim(unique(dat[,c(1,2)]))
# id1 is in increasing order
all(diff(dat[, 1]) >= 0)



# In total, 22646 unique sibling pairs. About 21338 unique ids
indx = which(r1 <= dat$IBS0 & dat$IBS0 <= r2 & dat$Kinship <= rk2 & dat$Kinship >= rk1)
sib = dat[indx, ]

# In total, 6263 rows
indx = which(dat$IBS0 <= r1 & dat$Kinship <= rk2 & dat$Kinship >= rk1)
trio = dat[indx, ]

# Creat parent offspring pair
family3 = c()
id_unique = unique(c(trio$ID1, trio$ID2))
for (i in 1:length(id_unique)){
  id = id_unique[i]
  sub = trio[which(trio$ID1 == id | trio$ID2 == id), ]
  if (nrow(sub) == 2){
    family_ids = unique(c(sub$ID1, sub$ID2))
    par_id = family_ids[family_ids != id]
    fam = c(par_id, id)
    family3 = rbind(family3, fam)
  }
}

# In total, get 1320 families. UKBiobank paper reports 1,066 sets of trios (two parents and an offspring)
# See https://www.nature.com/articles/s41586-018-0579-z
colnames(family3) = c("ID.par1", "ID.par2", "ID.child")
rownames(family3) = NULL



