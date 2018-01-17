#!/bin/bash

cd $REDUX_ENV

if [ "$1" == "clades" ]; then
    cd ../clades
    DB="clades"
else
    #cd redux2
    DB="variant"
fi

tbls=`sqlite3 $DB.db ".tables"`
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
    cnt=`sqlite3 $DB.db "select count(*) from $t"`
    echo "$t: $cnt"
done
echo "--------------------"
echo ""
