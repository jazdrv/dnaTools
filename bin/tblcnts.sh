#!/bin/bash

cd $REDUX_ENV

tbls=`sqlite3 variant.db ".tables"`
tbls1=`echo $tbls|tr " " "\n" | sort -g`
#echo "--------------------"
#echo "Tables"
#echo "--------------------"
#echo $tbls
#echo "--------------------"

echo ""
echo "--------------------"
echo "Table counts"
echo "--------------------"
for t in `echo "$tbls1"`; do
    cnt=`sqlite3 variant.db "select count(*) from $t"`
    echo "$t: $cnt"
done
echo "--------------------"
echo ""
