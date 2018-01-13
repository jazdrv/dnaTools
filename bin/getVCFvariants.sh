#!/bin/bash

#Note: bash hack for getting out variants and their quality data from the VCF files:
#This extracts all variants where more than 10% of the reads are called for that allele (accounting for 
#multiple possibilities for the variant at a given position).
#ie: ./getVCFvariants.sh "unzip/Donald-120386.vcf"

FILE="$1"
FILE="$REDUX_ENV/unzips/Newell-345238.vcf"
if [ "$1" == "" ]; then
    echo "Nothing given. Aborting."
    exit
fi
if [[ ! -f "$FILE" ]]; then
    echo "$FILE doesn't exist. Aborting."
    exit
fi

awk 'substr($0,1,1)!="#" {split($8,a,";"); split($10,b,":"); if ($7=="PASS") {c="PASS"} else {c="FAIL"}; print $2,$4,$5,c,a[1],a[6],b[2],b[3]}' \
"$FILE" |\
sed 's/=/ /g' |\
awk '$9~"," {sub(/^[0-9]*,/, "", $9)} {print $1,$2,$3,$4,$6,$8,$10,$9}' |\
awk '$0!~"," {$8/=$7;print} $0~"," {n=split($3,a,",");split($8,b,",");for(i=1;i<=n;i++) print $1,$2,a[i],$4,$5,$6,$7,b[i]/$7}' |\
awk '$8>0.1'

