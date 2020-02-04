# 4-2-2020 JHZ

library(gap)
cwd <- getwd()
source("swath-ms.ini")
load("swath-ms.rda")

# protein data

# EDA plots
df <- prot
rownames(df) <- gsub("DIINT01_","",rownames(df))
pdf("scatter-histogram-boxwhisker.pdf")
par(mfrow=c(3,1))
sapply(seq(1, ncol(df)), plotfun)
dev.off()

# Box-Whisker plot
pdf("box.pdf", width=50, height=12)
np <- ncol(df)
xtick <- seq(1, np, by=1)
boxplot(df,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(df), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("raw data")
sdf <- scale(df,scale=FALSE)
boxplot(sdf,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(sdf), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("mean-centred data")
ssdf <- scale(df,scale=TRUE)
boxplot(ssdf,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(ssdf), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("mean-centred and scaled data")
dev.off()

# outliers by AE
prot <- swath_protein[,-(1:2)]
r <- ae_swath(prot,hidden.layers=c(100,20,30))
idr <- cbind(swath_protein[c("Internal.ID","External.ID")],mse=rowSums(r)/ncol(prot))
ord <- with(idr, order(mse,decreasing=TRUE))
head(idr[ord,],30)
pdf("ae.pdf")
plot(idr[ord,3],cex=0.4)
dev.off()

# UMAP
library(uwot)
pdf("umap.pdf")
df_umap <- umap(prot,n_neighbors = 5, learning_rate = 0.5, pca=50, spread = 3)
plot(df_umap,col="black",cex=0.4)
outliers <- subset(prot,rownames(prot)%in%head(idr[ord,2],13))
points(df_umap[(1:196)[head(ord,13)],],col="blue",cex=0.4)
title(paste("UMAP with n_neighbors=5, learning_rate=0.5, pca=50, spread=3"))
dev.off()

# PCA
ppc <- with(prot, prcomp(na.omit(prot), rank=50, scale=TRUE))
write.table(data.frame(External.ID=rownames(prot),with(ppc,x)),file="pca.txt",sep="\t",quote=FALSE,row.names=FALSE)
pdf("pca.pdf")
screeplot(ppc, npcs=20, type="lines", main="PCA screeplot")
with(ppc, {
  plot(x[,1:2], main="PC1 -- PC2", cex=0.5, pch=21)
  legend("bottom", legend = c("lower detection", "lipemic"), box.lty = 0, cex = 0.8,
         col = c("red", "blue"), horiz = TRUE, inset = c(0, 1), xpd = TRUE, pch = c(10, 16))
})
biplot(ppc,cex=0.1)
title("biplot")
dev.off()

# RLE plot
pdf("rle.pdf", width=50, height=12)
par(mfrow=c(2,1))
makeRLEplot(prot, log2.data=FALSE, groups=group, col.group=col.group, cex=0.3, showTitle=TRUE)
makeRLEplot(prot, log2.data=FALSE, cex=0.3, showTitle=TRUE, title="Uncoloured relative log expression (RLE) plot")
dev.off()

# cluster analysis
km <- kmeans(prot,3)
with(km, table(cluster))
library(mclust)
mc <- Mclust(prot)
summary(mc)
pdf("mc.pdf")
plot(mc, what = "BIC")
dev.off()

# correlation
with(pheno_protein,cor(crp,CRP,use="complete.obs"))
with(pheno_protein,cor(transf,TRFE,use="complete.obs"))

# regression
attach(pheno_protein)
options(width=500)
cat("protein","intercept","sex",sep="\t",file="sex.tsv")
cat("\n",append=TRUE,file="sex.tsv")
cat("protein","intercept","sex","age","bmi",sep="\t",file="lm.tsv")
cat("\n",append=TRUE,file="lm.tsv")
sapply(seq(1, ncol(df)), regfun)
detach(pheno_protein)
pdf("qq.pdf")
sex <- read.delim("sex.tsv",as.is=TRUE)
gap::qqunif(with(sex,sex),cex=0.4)
title("QQ plot for sex")
lm <- read.delim("lm.tsv",as.is=TRUE)
gap::qqunif(with(lm,bmi),cex=0.4)
title("QQ plot for BMI")
dev.off()

Caprion <- paste(Sys.getenv("INF"),"Caprion",sep="/")
setwd(Caprion)
load("caprion.rda")

# genetics and phenotypes for association analysis

affymetrix <- function()
{
  affymetrix.id <- with(phenotypes,affymetrix_gwasqc_bl)
  write.table(affymetrix.id[!is.na(affymetrix.id)],file="affymetrix.id",row.names=FALSE,col.names=FALSE,quote=FALSE)
  uniprot <- scan("SomaLogic.uniprot",what="")
  SomaLogic <- subset(pap,Accession%in%uniprot)
  d1 <- t(SomaLogic[,-c(1:3)])
  d1 <- data.frame(caprion_id=rownames(d1),round(d1,3))
  protein <- merge(phenotypes[c("caprion_id","affymetrix_gwasqc_bl")],d1,by="caprion_id")
  p <- protein[,-1]
  p <- within(p,{for(i in names(p)[-1]) assign(paste0(i,"_invn"), gap::invnormal(p[i]))})
  p <- p[,-ncol(p)]
  id1_id2_missing_covariates_phenotypes <- merge(id1_id2_missing_covariates,p,
                                                 by.x="ID_1",by.y="affymetrix_gwasqc_bl")
  snptest_sample(id1_id2_missing_covariates_phenotypes,"SomaLogic.sample",
                 C=c("age","bmi",paste0("PC",1:20)),
                 D="sex",
                 P=names(p)[-1])
}

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

tromso_sample <- function()
# replication of TROMSO study
{
  d <- tromso_xlsx()
  pp <- names(table(with(d,sheet3["Protein.Name"])))
  overlap <- colnames(t1)[colnames(t1)%in%pp]
  selected <- with(d, sheet3[["Protein.Name"]] %in% overlap)
  s <- with(d, sheet3)[selected,c("Protein.Name","Ensembl.ID","Sentinel.Variant")]
  rsid <- levels(as.factor(s$Sentinel.Variant))
  ps <- phenoscanner::phenoscanner(snpquery=rsid, catalogue="pQTL")
  snps <- with(ps,snps)
  print(sort(as.numeric(levels(with(snps,chr)))))
  tromso <- merge(s,snps,by.x="Sentinel.Variant",by.y="snp")
  write.table(tromso[c("Sentinel.Variant","Protein.Name","chr")],
              file="tromso.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
  prot <- pheno_protein[c("caprion_id","affymetrix_gwasqc_bl",overlap)]
  prot <- within(prot,{for(i in names(prot[,-(1:2)])) assign(paste0(i,"_invn"),gap::invnormal(prot[i]))})
  prot <- prot[,-ncol(prot)]
  id1_id2_missing_covariates_phenotypes <- merge(id1_id2_missing_covariates,prot[,-1],
                                                 by.x="ID_1",by.y="affymetrix_gwasqc_bl")
  m <- merge(id1_id2_missing[,-3],id1_id2_missing_covariates_phenotypes,by=c("ID_1","ID_2"),all.x=TRUE)
  snptest_sample(m,"tromso.sample",
                 C=c("age","bmi",paste0("PC",1:20)),
                 D="sex",
                 P=names(prot[,-(1:2)]))
}

peptides_sample <- function()
# analysis on peptides
{
  p <- "ERAP2"
  d <- extract_peptide("ERAP2")
  id <- pheno_protein[,1:2]
  peptides <- merge(id,d,by="caprion_id")
  peptides <- within(peptides,{for(i in names(peptides[,-(1:2)])) assign(paste0(i,"_invn"),gap::invnormal(peptides[i]))})
  peptides <- peptides[,-ncol(peptides)]
  id1_id2_missing_covariates_phenotypes <- merge(id1_id2_missing_covariates,peptides[,-1],
                                                 by.x="ID_1",by.y="affymetrix_gwasqc_bl")
  snptest_sample(id1_id2_missing_covariates_phenotypes,paste0(p,".sample"),
                 C=c("age","bmi",paste0("PC",1:20)),
                 D="sex",
                 P=names(peptides)[-(1:2)])
}

