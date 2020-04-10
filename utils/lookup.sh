#!/usr/bin/bash

function check()
{
  export prot=$(grep -w $1 swath-ms.merge | cut -f5 | sed 's/_invn//g')
  if [ "${prot}" == "" ]; then
    echo Empty
  else
    grep -H -w ${prot} $INF/doc/hgTables.tsv
    grep -H -w INTERVAL-box.tsv
  fi
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

function Sun()
{
  awk 'NR>1{gsub(/_invn/,"");print $5,$6}' swath-ms.merge | \
  parallel -C' ' '
    echo {1} {2}
    grep -w {1} pQTL.Sun-B_pQTL_EUR_2017 | grep {2}
    grep -H -w {1} INTERVAL_box.tsv
  '
}

Sun
