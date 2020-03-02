# 2-3-2020

INF <- Sys.getenv("INF")

CRP <- within(read.table(paste(INF,"Caprion","plink2","CRP-plink2.gz",sep="/"),comment.char="",header=TRUE),{CRP=-log10(P)})
head(CRP)
P02741 <- within(read.table(paste(INF,"SWATH-MS","plink2","P02741-plink2.gz",sep="/"),comment.char="",header=TRUE),{P02741=-log10(P)})
head(P02741)

m <- merge(CRP,P02741,by=c("X.CHROM","POS"))
with(m,plot(CRP,P02741))
