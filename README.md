# SWATH-MS work

Acronym of Sequential Window Acquisition of all Theoretical mass spectrometry (SWATH-MS) / data independent acquisition (DIA)

* [swath-ms.R](swath-ms.R) ([swath-ms.ini](swath-ms.ini)) and [swath-ms.sh](swath-ms.sh) perform data proceessing and exploratory analysis. From this, we have
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
* [swath-ms.ipynb](swath-ms.ipynb) is a Jupyter notebook with some preprocessing done by [utils/tensorqtl.sh](utils/tensorqtl.sh).

Therefore with the eventual option of pilot samples plut genotypes with an MAF cutoff 0.01, SLURM may or may not be needed.

## References

Ludwig C, Gillet L, Rosenberger G, Amon S, Collins BC, Aebersold R (2018). Data‐independent acquisition‐based SWATH‐MS for quantitative proteomics: 
a tutorial. Mol Syst Biol 14:e8126, https://doi.org/10.15252/msb.20178126.

[OpenSWATH](http://openswath.org/en/latest/),
[ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz).
