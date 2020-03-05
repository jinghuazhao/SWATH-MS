# SWATH-MS work

* [swath-ms.R](swath-ms.R) ([swath-ms.ini](swath-ms.ini)) performs data proceessing and exploratory analysis. From this, we have
  * [utils/qctool.sh](utils/qctool.sh) attempts to extract genotype/sample information into a single .bgen file gives way to [qctool.sb](utils/qctool.sb) for its speed.
  * [utils/quicktest.sb](utils/quicktest.sb) is SLURM version of QUICKTEST.
  * [utils/quicktest.sh](utils/quicktest.sh) is non-SLURM version of QUICKTEST.
  * [utils/plink2.sb](utils/plink2.sb) is SLURM version of PLINK2.
  * [utils/plink2.sh](utils/plink2.sh) non-SLURM version of PLINK2.
  * [utils/check.R](utils/check.R) attempts to check for consistence with Caprion.
  * [utils/affymetrix.sh](utils/affymetrix.sh) performs the association analysis on specific protein-variant combinations.
  * [utils/pgwas.sh](utils/pgwas.sh) ([utils/pgwas.ini](utils/pgwas.ini)) conducts genomewide associations analyses on all proteins. Note that
    * BOLT-LMM took 24hr for one protein from all data (N=43,059) but failed to run on available genotypes and samples (N=196). It also uses 8-bit version of bgen (qctool -bgen-bits 8 and also the master genotype files).
    * PLINK2 is attractive with its speed.
    * QUICKTEST is faster than SNPTEST and takes into account uncertainty, and is therefore more preferable.
    * SNPTEST gives verbose screen output with -printids.
  * [utils/sentinels_nold.sh](utils/sentinels_nold.sh) and [utils/merge.sh](utils/merge.sh) performs sentinel selection.
* [swath-ms.ipynb](swath-ms.ipynb) is a Jupyter notebook calling tensorQTL.

Therefore with the eventual option of pilot samples plut genotypes with an MAF cutoff 0.01, SLURM may or may not be needed.

## URLs

[SWATH-MS](https://imsb.ethz.ch/research/aebersold/research/swath-ms.html),
[OpenSWATH](http://openswath.org/en/latest/),
[ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz).
