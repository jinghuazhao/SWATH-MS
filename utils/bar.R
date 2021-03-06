# 11-12-2019 JHZ

bar <- function(data, title)
{
  affy <- read.delim(data,as.is=TRUE)
  affy_vars <- c("chromosome","position","add_beta_1","add_se_1","add_pvalue","A","B","uniprot")
  box <- read.table("SomaLogic.box",as.is=TRUE,header=TRUE,fill=TRUE)
  box <- within(box,{Allele1=toupper(Allele1);Allele2=toupper(Allele2)})
  box_vars <- c("chr","pos","Effect","StdErr","logP","protein","SOMAS_ID_round1","rsid","Allele1","Allele2")
  affybox <- merge(box[box_vars],affy[affy_vars],by.x=c("chr","pos"),by.y=c("chromosome","position"))
  affybox <- within(affybox[-5,],{switch <- B!=Allele1;add_beta_1[switch] <- -add_beta_1[switch]})
  SomaLogic <- data.frame(affybox[c("Effect","StdErr","logP","protein","uniprot","rsid")],platform="SomaLogic")
  swath <- data.frame(affybox[c("add_beta_1","add_se_1","add_pvalue","protein","uniprot","rsid")],platform="SWATH")
  names(swath)[1:3] <- c("Effect","StdErr","logP")
  swath <- within(swath, {logP <- log10(logP)})
  tabbedData <- within(rbind(SomaLogic,swath),{protein <- paste(protein,rsid,sep="-")})
  limits <- with(tabbedData, aes(ymax = Effect + StdErr, ymin = Effect - StdErr))
  p <- ggplot(data = tabbedData, aes(x = factor(protein), y = Effect, fill = factor(platform)))
  p + geom_bar(stat = "identity", position = position_dodge(0.9)) + geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) + labs(x = "Protein-SNP pair", y = "Effect size") + ggtitle(title) + scale_fill_discrete(name = "platform") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

library(ggplot2)
pdf("bar.pdf")
bar("affymetrix.tsv","Effect size by protein and platform (SWATH untransformed)")
bar("affymetrix_invn.tsv","Effect size by protein and platform (SWATH invnorm)")
dev.off()
