# 28-2-2020 JHZ

export interval=/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
qctool -g $interval/impute_#_interval.bgen -s $interval/interval.samples -incl-samples affymetrix.id \
       -og swath-ms.bgen -os swath-ms.samples

qctool -g swath-ms-#.bgen -s swath-ms.sample -ofiletype binary_ped -og swath-ms.bgen
