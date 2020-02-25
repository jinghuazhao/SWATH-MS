# SWATH-MS work

* [qctool.sh](utils/qctool.sh) attempts to extract genotype/sample information into a single .bgen file but is too slow.
* [swath-ms.ini](swath-ms.ini) and [swath-ms.R](swath-ms.R) perform data proceessing.
* [affymetrix.sh](utils/affymetrix.sh) performs the association analysis on specific protein-variant combinations.
* [pgwas.sh](utils/pgwas.sh) conducts genomewide associations analyses on all proteins.
  * BOLT-LMM took 24hr for one protein from all data (N=43,059) but failed to run on available genotypes and samples (N=196). It also uses 8-bit version of bgen (qctool -bgen-bits 8 and also the master genotype files).
  * SNPTEST gives verbose screen output with -printids.
  * QUICKTEST is faster than SNPTEST.

## URLs

[SWATH-MS](https://imsb.ethz.ch/research/aebersold/research/swath-ms.html),
[OpenSWATH](http://openswath.org/en/latest/),
[ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz).
