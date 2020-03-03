# 3-3-2020 JHZ

export interval=/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
# SWATH-MS samples
qctool -g $interval/impute_#_interval.bgen -s $interval/interval.samples -incl-samples affymetrix.id \
       -og swath-ms.bgen -os swath-ms.samples
# MAF cutoff 0.01
qctool -g swath-ms-#.bgen -s swath-ms.sample -ofiletype binary_ped -og swath-ms.bgen
plink --bfile swath-ms.bgen --maf 0.01 --make-bed --out swath-ms.01
awk '{$1=$2};1' swath-ms.01.fam > swath-ms.fam
cut -f2 swath-ms.01.bim > swath-ms.01.snpids
qctool -g swath-ms-#.bgen -s SomaLogic.sample -og swath-ms.01.bgen -bgen-bits 8 -incl-snpids swath-ms.01.snpids
bgenix -g swath-ms.01.bgen -index -clobber

