#!/usr/bin/bash

export TMPDIR=/rds/user/$USER/hpc-work/

module load plink/2.00-alpha

function binary_ped()
{
seq 609 | \
parallel -C' ' '
  export col=$(cut -d" " -f {} swath-ms.uniprot); \
  for v in ${col} ${col}_invn;
  do
      echo ${v}
      plink2 \
             --bfile swath-ms.01 --fam swath-ms.fam \
             --glm hide-covar --input-missing-phenotype -9 --covar-variance-standardize \
             --pheno swath-ms.pheno --pheno-name ${v} --covar swath-ms.covar \
             --out work/${v}
      grep -v NA work/${v}.${v}.glm.linear | \
      gzip -f > plink2/${v}-plink2.gz
      rm work/${v}.${v}.glm.linear
  done
'
}

function bgen()
{
seq 609 | \
parallel -C' ' '
  export col=$(cut -d" " -f {} swath-ms.uniprot); \
  for v in ${col} ${col}_invn;
  do
      echo ${v}
      plink2 \
             --bgen swath-ms.01.bgen --sample swath-ms.sample \
             --glm hide-covar --input-missing-phenotype -9 --covar-variance-standardize \
             --pheno swath-ms.pheno --pheno-name ${v} --covar swath-ms.covar \
             --out work/${v}
      grep -v NA work/${v}.${v}.glm.linear | \
      gzip -f > bgen/${v}-plink2.gz
      rm work/${v}.${v}.glm.linear
  done
'
}

bgen
