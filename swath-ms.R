# 3-2-2020 JHZ

cwd <- getwd()
library(gap)
source("swath-ms.ini")
load("swath-ms.rda")

Caprion <- paste(Sys.getenv("INF"),"Caprion",sep="/")
setwd(Caprion)
load("caprion.rda")
names(Samples) <- c("caprion_id","external_id","comment")
phenotypes <- read.delim("interval_caprion_pilot_samples_phenotype_data.tsv",as.is=TRUE)
phenotypes <- within(phenotypes,{
  age <- agepulse
  sex <- sexpulse
  bmi <- wt_bl/ht_bl/ht_bl
  crp <- crp_bl
  transf <- transf_bl
})
pd <- merge(phenotypes[c("caprion_id","affymetrix_gwasqc_bl","sex","age","bmi","crp","transf")],Samples,by="caprion_id")

id1_id2_0 <- read.table("interval.samples",skip=2,col.names=c("ID_1","ID_2","missing"))
missing <- read.table("merged_imputation.missing",col.names=c("affymetrix_gwasqc_bl","missing"))
id1_id2_missing <- merge(id1_id2_0[,-3],missing,by.x="ID_1",by.y="affymetrix_gwasqc_bl")
eigenvec <- read.delim("merged_imputation.eigenvec")
covariates <- merge(pheno_protein[c("affymetrix_gwasqc_bl","sex","age","bmi")],eigenvec[,-1],
                    by.x="affymetrix_gwasqc_bl",by.y="IID")
id1_id2_missing_covariates <- merge(id1_id2_missing,covariates,by.x="ID_1",by.y="affymetrix_gwasqc_bl")
setwd(cwd)
