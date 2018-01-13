#!/bin/bash

#Note: temporary script ... used with getVCFVariants.sh to push that data into sqlite3 

cd $REDUX_ENV

echo "create table variants_sh (variant int,tot int,pass tinyint,fail tinyint);" > $REDUX_ENV/go.out
cat $REDUX_ENV/data.out|awk 'BEGIN{p=0;f=0;c=0;k="beg"}{
  if (k != $1) {
    if (k != "beg") {
        print "insert into variants_sh values("k","c","p","f");";
    }
    k=$1;p=0;f=0;c=0;
  }
  if ($0 ~ /PASS/) {p++;c++;};
  if ($0 ~ /FAIL/) {f++;c++;};
}' >> $REDUX_ENV/go.out

