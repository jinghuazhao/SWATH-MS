#!/usr/bin/bash

export TMPDIR=/rds/user/$USER/hpc-work/

seq 609 | \
parallel -C' ' '
  export col=$(cut -d" " -f {} swath-ms.uniprot)
  quicktest --geno swath-ms.01.bgen --pheno swath-ms.sample --snptest --bgen \
            --npheno ${col} --ncovar sex --ncovar age --ncovar bmi \
            --ncovar PC1 --ncovar PC2 --ncovar PC3 --ncovar PC4 --ncovar PC5  \
            --ncovar PC6 --ncovar PC7 --ncovar PC8 --ncovar PC9 --ncovar PC10 \
            --ncovar PC11 --ncovar PC12 --ncovar PC13 --ncovar PC14 --ncovar PC15  \
            --ncovar PC16 --ncovar PC17 --ncovar PC18 --ncovar PC19 --ncovar PC20 \
            --method-mean --method-score --compute-rSqHat --out work/${col}-qt
  grep -v NA work/${col}-qt | \
  gzip -f > quicktest/${col}-qt.gz
  rm work/${col}-qt
'

# _invn version:
# export col=$(cut -d" " -f {} swath-ms.uniprot)_invn
