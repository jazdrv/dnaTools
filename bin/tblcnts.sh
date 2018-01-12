#!/bin/bash

cd $REDUX_ENV

tbls=`
sqlite3 variant.db<<!
.tables
!`
echo $tbls

for t in `echo $tbls`; do
sqlite3 variant.db<<!
select count(*) from $t;
!
done
