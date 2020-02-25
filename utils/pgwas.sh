# 25-2-2020 JHZ

source utils/pgwas.ini

# initialize

for i in `seq 1 3`; do pgwas_quicktest; done

(
  cat *out | head -19 | sed 's/allele//g;s/frequentist_//g' | tail -n 1 | awk -v OFS="\t" '{print "uniprot", "protein", $0}'
  for i in `seq 0 3`
  do
    awk 'NR==20' ${uniprot[i]}-${rsid[i]}.out | \
    awk -v uniprot=${uniprot[i]} -v protein=${protein[i]} -v OFS="\t" '{print uniprot, protein, $0}'
  done
) > affymetrix.tsv
(
  cat *out | head -19 | sed 's/allele//g;s/frequentist_//g' | tail -n 1 | awk -v OFS="\t" '{print "uniprot", "protein", $0}'
  for i in `seq 0 3`
  do
    awk 'NR==20' ${uniprot[i]}_invn-${rsid[i]}.out | \
    awk -v uniprot=${uniprot[i]} -v protein=${protein[i]}_invn -v OFS="\t" '{print uniprot, protein, $0}'
  done
) > affymetrix_invn.tsv

R --no-save -q < $Caprion/utils/ps.R

(
  head -1 $Caprion/SOMALOGIC_Master_Table_160410_1129info.tsv
  cut -f1 affymetrix.tsv | sed '1d' | grep -w -f - $Caprion/SOMALOGIC_Master_Table_160410_1129info.tsv
) > SomaLogic.info

join -j2 <(cut -f2,10 SomaLogic.info | sed '1d'| sort | uniq) \
         <(cut -f10 SomaLogic.info | sed '1d'| sort | uniq | grep -f - -w $Caprion/glist-hg19 | cut -d' ' -f1,4) | \
cut -d' ' -f2,3 | \
parallel -j1 -C' ' 'cd $Caprion/SomaLogic;unzip -o {1}.zip {1}/{1}_chrom_{2}_meta_final_v1.tsv.gz;cd -'

(
  echo "SOMAS_ID_round1 chr pos rsid protein VARIANT_ID chromosome position Allele1 Allele2 Effect StdErr logP"
  join -j2 <(cut -f2,10 SomaLogic.info | sed '1d'| sort | uniq) \
           <(cut -f10 SomaLogic.info | sed '1d'| sort | uniq | grep -f - -w $Caprion/glist-hg19 | cut -d' ' -f1,4) | \
  join - <(cut -f2,3,6 affymetrix.tsv | sed '1d' | sort -k1,1) | cut -d' ' -f1-5 | \
  parallel -j1 -C' ' 'echo {2} {3} {5} {4} {1} $(zgrep -w {5} $Caprion/SomaLogic/{2}/{2}_chrom_{3}_meta_final_v1.tsv.gz)'
) > SomaLogic.box

R --no-save -q <<END
  s <- read.delim("SomaLogic.info",as.is=T)
  v <- c("chromosome_number","start_position","end_position","ensembl_gene_id","external_gene_name","UniProt","Target")
  s <- within(s[v],{chromosome_number <- paste0("chr",chromosome_number)})
  names(s) <- c("chrom","chromStart","chromEnd","ENSEMBL","HGNC","uniprot","Protein")
  write.table(s,file="SomaLogic.gene",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
END

(
  echo $(head -1 SomaLogic.gene | sed 's/chrom/geneChr/;s/chromStart/geneStart/;s/chromEnd/geneEnd/') \
       $(head -1 affymetrix.results | sed 's/chrom/snpChr/;s/chromStart/snpStart/;s/chromEnd/snpEnd/') | \
  sed 's/ /\t/g'
  bedtools intersect -a SomaLogic.gene -b affymetrix.results -wo
) | \
cut -f1-20 > SomaLogic.tsv

(
  echo chrom chromStart chromEnd gene
  cut -f2 affymetrix.tsv | \
  sed '1d' | \
  grep -w -f - $Caprion/glist-hg19 | \
  awk -v OFS="\t" '{print "chr" $1,$2,$3,$4}'
) | \
bedtools intersect -a - -b affymetrix.results -wo 

sed 's/Caprion/SWATH/g;s/caprion/swath/g'  $Caprion/utils/bar.R > utils/bar.R
R --no-save -q < utils/bar.R
