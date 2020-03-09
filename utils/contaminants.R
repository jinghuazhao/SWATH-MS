# 9-3-2020 JHZ

source("swath-ms.ini")

# original Uniprot IDs
u <- read.table("contaminants.uniprot",as.is=TRUE,col.name="Uniprot")
library(openxlsx)
rds.xlsx <- "DIINT01_BATCH_1_MS_BATCH_INPUT_RDS_02DEC19 final.xlsx"
protein_matrix.file <- "DIINT_protein_matrix final.xlsx"
protein_matrix <- read.xlsx(protein_matrix.file, sheet = 1, startRow = 1)
protein_contaminants <- merge(protein_matrix,u,by.x="Internal.ID",by.y="Uniprot")
rownames(protein_contaminants) <- protein_contaminants[["Internal.ID"]]
with(protein_contaminants,Internal.ID)

# EDA plots
df <- t(protein_contaminants[,-1])
pdf("scatter-histogram-boxwhisker-contaminants.pdf")
par(mfrow=c(3,1))
sapply(seq(1, ncol(df)), plotfun)
dev.off()

## EBSEMBL + RefSeq lookup
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
attr <- listAttributes(ensembl)

# ENSEMBL
b <- c("btaurus_homolog_ensembl_peptide","btaurus_homolog_canonical_transcript_protein","btaurus_homolog_ensembl_gene")
e <- c("ensembl_peptide_id","ensembl_gene_id")
bm_be <- getBM(attributes = c(b,e), mart=ensembl)
es <- scan("contaminants.ensembl","")
s <- subset(bm_be,bm_be[["btaurus_homolog_ensembl_peptide"]]%in%es|bm_be[["btaurus_homolog_canonical_transcript_protein"]]%in%es)
unique(s[["btaurus_homolog_ensembl_peptide"]])
write.table(s,file="contaminants_annot.txt",row.names=FALSE)

# RefSeq
rs <- scan("contaminants.refseq","")
r <- c("refseq_peptide","refseq_mrna","refseq_peptide_predicted")
bm_r <- getBM(attributes=r, mart=ensembl)
subset(bm_r,bm_r[["refseq_peptide_predicted"]]%in%rs)
