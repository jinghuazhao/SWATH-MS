# SWATH-MS work

* [swath-ms.ini](swath-ms.ini) and [swath-ms.R](swath-ms.R) perform data proceessing.
  * [utils/qctool.sh](utils/qctool.sh) attempts to extract genotype/sample information into a single .bgen file gives way to [qctool.sb](utils/qctool.sb) for its speed.
  * [utils/affymetrix.sh](utils/affymetrix.sh) performs the association analysis on specific protein-variant combinations.
  * [utils/pgwas.sh](utils/pgwas.sh) conducts genomewide associations analyses on all proteins. Note that
    * BOLT-LMM took 24hr for one protein from all data (N=43,059) but failed to run on available genotypes and samples (N=196). It also uses 8-bit version of bgen (qctool -bgen-bits 8 and also the master genotype files).
    * PLINK2 is reasonable.
    * SNPTEST gives verbose screen output with -printids.
    * QUICKTEST is faster than SNPTEST and takes into account uncertainty.

QUICKTEST is therefore the chosen software.

## URLs

[SWATH-MS](https://imsb.ethz.ch/research/aebersold/research/swath-ms.html),
[OpenSWATH](http://openswath.org/en/latest/),
[ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz).
