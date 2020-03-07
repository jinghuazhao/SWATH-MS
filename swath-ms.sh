# 7-3-2020 JHZ

grep \> contaminants.fasta  | cut -d' ' -f1 | awk 'length($1)<=8'| sed 's/>//g' > contaminants.uniprot
