# 18-2-2020 JHZ

cwd <- getwd()
Caprion <- paste(Sys.getenv("INF"),"Caprion",sep="/")
setwd(Caprion)
source("caprion.ini")
setwd(cwd)
source("swath-ms.ini")
swath_ms_data()
load("swath-ms.rda")
swath_overlap()
load("swath-ms-overlap.rda")
prot <- swath_protein[,-(1:2)]

# Outliers by AE
# module load python/3.6
# source $HOME/.conda/tensorflow-gpu/bin/activate
# protein data

mse_threshold <- 0.0006
r <- ae_swath(prot,hidden.layers=c(800,150,400))
idr <- cbind(swath_protein[c("Internal.ID","External.ID")],mse=rowSums(r)/ncol(prot))
ord <- with(idr, order(mse,decreasing=TRUE))
print(subset(idr[ord,],mse>mse_threshold),row.names=FALSE)
excl <- with(subset(idr,mse>mse_threshold),External.ID)
outliers <- rownames(prot)%in%excl
pdf("ae.pdf")
plot(idr[ord,3],cex=0.4)
dev.off()

library(gap)

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

# UMAP
library(uwot)
pdf("umap.pdf")
df_umap <- umap(prot,n_neighbors = 5, learning_rate = 0.5, pca=50, spread = 3)
plot(df_umap,col="black",cex=0.4)
points(df_umap[outliers,],col="red",cex=0.4)
title(paste("UMAP with n_neighbors=5, learning_rate=0.5, pca=50, spread=3"))
dev.off()

# PCA
ppc <- with(prot, prcomp(na.omit(prot), rank=50, scale=TRUE))
write.table(data.frame(External.ID=rownames(prot),with(ppc,x)),file="pca.csv",sep=",",quote=FALSE,row.names=FALSE)
pdf("pca.pdf")
screeplot(ppc, npcs=20, type="lines", main="PCA screeplot")
with(ppc, {
  plot(x[,1:2], main="PC1 -- PC2", cex=0.5, pch=21)
  points(x[outliers,1:2],cex=0.5,pch=21,col="red")
  legend("bottom", legend = c(paste("rms >", mse_threshold), "other"), box.lty = 0, cex = 0.8,
         col = c("red", "black"), horiz = TRUE, inset = c(0, 1), xpd = TRUE, pch = c(10, 16))
})
biplot(ppc,cex=0.1)
title("biplot")
dev.off()

# RLE plot
pdf("rle.pdf", width=50, height=12)
par(mfrow=c(2,1))
groups <- rep(1,nrow(prot))
groups[outliers] <- 2
col.groups=c("black","red")
makeRLEplot(t(prot), log2.data=FALSE, groups=groups, col.group=col.groups, cex=0.3, showTitle=TRUE)
makeRLEplot(t(prot), log2.data=FALSE, cex=0.3, showTitle=TRUE, title="Uncoloured relative log expression (RLE) plot")
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
CRP_TRFE <- swath_protein[c("Internal.ID","P02741","P02787")]
names(CRP_TRFE)[2:3] <- c("CRP","TRFE")
j <- merge(swath_pheno,CRP_TRFE,by.x="swathMS_id",by.y="Internal.ID")
with(j,cor(CRP_bl,CRP,use="complete.obs"))
with(j,cor(TRANSF_bl,TRFE,use="complete.obs"))
jj <- merge(j,pheno_protein,by="caprion_id")
with(jj,cor(CRP.x,CRP.y))
with(jj,cor(TRFE.x,TRFE.y))

# regression
pheno_protein <- merge(swath_pheno,swath_protein,by.x="swathMS_id",by.y="Internal.ID")
df <- pheno_protein[,-(1:16)]
options(width=500)
cat("protein","intercept","sex",sep="\t",file="sex.tsv")
cat("\n",append=TRUE,file="sex.tsv")
cat("protein","intercept","sex","age","bmi",sep="\t",file="lm.tsv")
cat("\n",append=TRUE,file="lm.tsv")
sapply(seq(1, ncol(df)), regfun)
pdf("qq.pdf")
sex <- read.delim("sex.tsv",as.is=TRUE)
gap::qqunif(with(sex,sex),cex=0.4)
title("QQ plot for sex")
lm <- read.delim("lm.tsv",as.is=TRUE)
gap::qqunif(with(lm,bmi),cex=0.4)
title("QQ plot for BMI")
dev.off()

# genetics and phenotypes for association analysis

id1_id2_0 <- read.table(paste(Caprion,"interval.samples",sep="/"),skip=2,col.names=c("ID_1","ID_2","missing"))
missing <- read.table(paste(Caprion,"merged_imputation.missing",sep="/"),col.names=c("affymetrix_gwasqc_bl","missing"))
id1_id2_missing <- merge(id1_id2_0[,-3],missing,by.x="ID_1",by.y="affymetrix_gwasqc_bl")
eigenvec <- read.delim(paste(Caprion,"merged_imputation.eigenvec",sep="/"))
covariates <- merge(pheno_protein[c("Affymetrix_gwasQC_bl","sex","age","bmi")],eigenvec[,-1],
                    by.x="Affymetrix_gwasQC_bl",by.y="IID")
id1_id2_missing_covariates <- merge(id1_id2_missing,covariates,by.x="ID_1",by.y="Affymetrix_gwasQC_bl")
affymetrix()
