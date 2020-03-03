# 3-3-2020 JHZ

export tag=_nold

function pgz()
# 1. extract all significant SNPs
{
  ls plink2/*.gz | grep -v inv| \
  sed 's|plink2/||g;s/.gz//g' | \
  parallel -j3 -C' ' '
  (
  # zcat plink2/{}.gz | head -1
    zcat plink2/{}.gz | awk "
    function abs(x)
    {
      if (x<0) return -x;
      else return x;
    }
    NR>1 && abs(\$11)>=6.496698" | sort -k1,1n -k2,2n
  ) | gzip -f > sentinels/{}.p.gz'
}

function _HLA()
# 2. handling HLA
{
  for p in $(ls plink2/*.gz | sed 's|plink2/||g;s/.gz//g')
  do
    (
      zcat plink2/${p}.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      zcat sentinels/${p}.p.gz | \
      awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
      awk '!($1 == "chr6" && $3 >= 25392021 && $3 < 33392022)'
      zcat sentinels/${p}.p.gz | \
      awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
      awk '$1 == "chr6" && $3 >= 25392021 && $3 < 33392022' | \
      sort -k13,13g | \
      awk 'NR==1'
    ) > sentinels/${p}${tag}.p
    export lines=$(wc -l sentinels/${p}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}${tag} with $lines lines
      rm sentinels/${p}${tag}.p
    fi
  done
}

for cmd in pgz _HLA; do $cmd; done
