#!/bin/bash

cd $REDUX_ENV

DB="variant"
info="c_variants|pos v1_hg19|grch37 v1_hg38|grch38 variants|pos"

for I in `echo $info`; do
tbl=`echo $I|cut -d'|' -f1`
fld=`echo $I|cut -d'|' -f2`
echo "-------------------"
echo $tbl
echo "-------------------"
echo "select * from $tbl order by $fld limit 5;"
sqlite3 $DB.db <<!
.header on
.mode column
select count(*) from $tbl;
select * from $tbl order by $fld limit 5;
!
done

for I in `echo $info`; do
tbl=`echo $I|cut -d'|' -f1`
fld=`echo $I|cut -d'|' -f2`
echo "-------------------"
echo $tbl
echo "-------------------"
echo "select * from $tbl order by $fld limit 5;"
sqlite3 $DB.db <<!
.schema $tbl
!
done

