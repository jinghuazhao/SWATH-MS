# 4-3-2020 JHZ

export uniprot=P08195
export protein=4F2
gunzip -c bgen/${uniprot}-plink2.gz | cut -f1,2,12 --output-delimiter=' ' > ${protein}.txt

R --slave --vanilla --args \
  input_data_path=${protein}.txt \
  output_data_rootname=${protein}_qq \
  plot_title="${protein} example" < turboqq.r

R --slave --vanilla --args \
  input_data_path=${protein}.txt \
  output_data_rootname=${protein}_man \
  reference_file_path=turboman_hg19_reference_data.rda \
  pvalue_sign=8.210181e-12 \
  plot_title="${protein} example" < turboman.r


# P04004 -  VTN
