# 9-3-2020 JHZ

grep \> contaminants.fasta  | cut -d' ' -f1 | awk 'length($1)<=8'| sed 's/>//g' > contaminants.uniprot
grep \> contaminants.fasta  | cut -d' ' -f1 | grep ENSEMBL | sed 's/>//g;s/ENSEMBL://g' > contaminants.ensembl
grep \> contaminants.fasta  | cut -d' ' -f1 | grep REFSEQ | sed 's/>//g;s/REFSEQ://g' > contaminants.refseq

R --no-save -q < utils/contaminants.R
R --no-save -q < swath-ms.R
