# 4-3-2020 JHZ

cat qqman.list | \
parallel -C' ' '
  echo {1} {2}
  export uniprot={1}
  export protein={2}
  gunzip -c bgen/${uniprot}-plink2.gz | \
  cut -f1,2,12 --output-delimiter=" " > ${uniprot}.txt

  R --slave --vanilla --args \
    input_data_path=${uniprot}.txt \
    output_data_rootname=${uniprot}_qq \
    plot_title="${uniprot[$i]} (${protein})" < turboqq.r

  R --slave --vanilla --args \
    input_data_path=${uniprot}.txt \
    output_data_rootname=${uniprot}_man \
    custom_peak_annotation_file_path=annotate.txt \
    reference_file_path=turboman_hg19_reference_data.rda \
    pvalue_sign=8.210181e-12 \
    plot_title="${uniprot} (${protein})" < turboman.r

  rm ${uniprot}.txt
'
