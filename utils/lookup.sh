#!/usr/bin/bash

function check()
{
  grep -H -w $1 swath-ms.merge
  export prot=$(grep $1 swath-ms.merge | cut -f5 | sed 's/_invn//g')
  grep -H -w ${prot} $INF/doc/hgTables.tsv
  grep -H -w ${prot} ~/SomaLogic/doc/SOMALOGIC_Master_Table_160410_1129info.tsv
}

function Sun()
{
  grep $1 -H -w pQTL.Sun-B_pQTL_EUR_2017
  check $1
}

function Folkersen()
{
  grep $1 -H -w pQTL.Folkersen-L_Proteins_EUR_2017
  check $1
}

function Suhre()
{
  grep $1 -H -w pQTL.Suhre-K_pQTL_EUR_2017
  check $1
}

function pQTL()
{
  grep $1 -H -w pQTL.pQTL_2017
  check $1
}

for rsid in $(sed '1d' pQTL.Sun-B_pQTL_EUR_2017 | awk '{print $2}' | uniq )
do
  echo --- ${rsid} ---
  Sun ${rsid}
  echo 
done
