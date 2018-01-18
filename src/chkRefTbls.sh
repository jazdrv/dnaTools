#!/bin/bash

cd $REDUX_ENV

DB="variant"

sqlite3 $DB.db <<!
.header on
.mode column
select * from c_variants order by pos limit 5;
select * from v1_hg19 order by grch37 limit 5;
select * from v1_hg38 order by grch38 limit 5;
!
