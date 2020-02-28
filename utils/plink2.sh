#!/usr/bin/bash

export TMPDIR=/rds/user/$USER/hpc-work/

module load plink/2.00-alpha

awk '{$1=$2};1' swath-ms.01.fam > swath-ms.fam
seq 609 | \
parallel -C' ' '
  export col=$(cut -d" " -f {} swath-ms.uniprot); \
  for v in ${col} ${col}_inv;
  do
      echo ${v}
      plink2 \
             --bfile swath-ms.01 --fam swath-ms.fam \
             --glm hide-covar --input-missing-phenotype -9 --covar-variance-standardize \
             --pheno swath-ms.pheno --pheno-name ${col} --covar swath-ms.covar \
             --out work/${v}
      grep -v NA work/${v}.${col}.glm.linear | \
      gzip -f > plink2/${v}-plink2.gz
      rm work/${v}.${col}.glm.linear
  done
'
