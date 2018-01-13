#!/bin/bash

cd $REDUX_ENV

tbls=`sqlite3 variant.db ".tables"`

#echo "--------------------"
#echo "Tables"
#echo "--------------------"
#echo $tbls
#echo "--------------------"

echo ""
echo "--------------------"
echo "Table counts"
echo "--------------------"
for t in `echo "$tbls"`; do
    cnt=`sqlite3 variant.db "select count(*) from $t"`
    echo "$t: $cnt"
done
echo "--------------------"
echo ""
