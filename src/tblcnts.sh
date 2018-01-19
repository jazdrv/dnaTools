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
tbls1=`echo $tbls|tr " " "\n" | sort -u`
#echo "--------------------"
#echo "Tables"
#echo "--------------------"
#echo $tbls
#echo "--------------------"

c_tbls1=`echo $tbls|tr " " "\n" | sort -u|grep "^c_"`
s_tbls1=`echo $tbls|tr " " "\n" | sort -u|grep "^s_"`
v1_tbls1=`echo $tbls|tr " " "\n" | sort -u|grep "^v1_"`
v2_tbls1=`echo $tbls|tr " " "\n" | sort -u|sed -e '/^c_/d' -e '/^s_/d' -e '/^v1_/d'`

echo ""
echo "--------------------"
echo "Table counts"
echo "--------------------"
for t in `echo "$v1_tbls1"`; do
    cnt=`sqlite3 $DB.db "select count(*) from $t"`
    echo "$t: $cnt"
done
for t in `echo "$c_tbls1"`; do
    cnt=`sqlite3 $DB.db "select count(*) from $t"`
    echo "$t: $cnt"
done
for t in `echo "$s_tbls1"`; do
    cnt=`sqlite3 $DB.db "select count(*) from $t"`
    echo "$t: $cnt"
done
for t in `echo "$v2_tbls1"`; do
    cnt=`sqlite3 $DB.db "select count(*) from $t"`
    echo "$t: $cnt"
done
echo "--------------------"
echo ""
