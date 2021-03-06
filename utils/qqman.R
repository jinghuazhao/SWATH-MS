# 2-4-2020 JHZ

uniprot <- Sys.getenv("uniprot");
print(uniprot);
gz <- gzfile(paste0("bgen/",uniprot,"-plink2.gz"));
require(qqman);
tbl <- read.delim(gz,as.is=TRUE);
tbl <- within(tbl,{
   SNP <- ID
   CHR <- as.numeric(X.CHROM)
   BP <- POS
})
tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
# qq <- paste0(uniprot,"_qq.png");
# png(qq,width=12,height=10,units="in",pointsize=4,res=300)
# qq(with(tbl,P))
# dev.off()
manhattan <- paste0(uniprot,"_manhattan.png");
png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
manhattan(tbl,main=uniprot,genomewideline=-log10(8.210181e-12),suggestiveline=FALSE,ylim=c(0,25));
dev.off();

# more check on 2/4/2020
# p04195 <- read.delim("P04196_invn-plink2.gz",as.is=TRUE)
# with(p04195,qq(P))
# d <- p04195[,c(1,2,12)]
# names(d) <- c("CHR","BP","P")
# manhattan(d)

