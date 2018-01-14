#!/bin/bash

#Note: temporary script ... used with getVCFVariants.sh to push that data into sqlite3 

cd $REDUX_ENV

#get the default hg38 data (doesn't need to be run every time)
#tail -n +2 $REDUX_ENV/snps_hg38.csv|cut -d'"' -f8,22,24 | sed -e 's/"/,/g' |sort -u> $REDUX_ENV/hg38.0

#echo "drop table variants_sh;" > $REDUX_ENV/go.tmp
#echo "create table variants_sh (allinfo char(20), code1 char(10), code2 char(10), variant int, tot tinyint, pass tinyint, fail tinyint);" \
#    >> $REDUX_ENV/go.tmp
cat $REDUX_ENV/new.0|awk 'BEGIN{p=0;f=0;t=0;a="beg";k="";c1="";c2="";}{
  if (a != $1"|"$2"|"$3) {
    if (a != "beg") {
        if (1==2) {
            print "insert into variants_sh values(+"a"+,+"c1"+,+"c2"+,"k","t","p","f");";
        }
        print a;
    }
    k=$1;c1=$2;c2=$3;
    a=k"|"c1"|"c2;
    p=0;f=0;c=0;
  }
  if ($0 ~ /PASS/) {p++;t++;};
  if ($0 ~ /FAIL/) {f++;t++;};
}' | sed -e "s/\+/'/g" |sort -u\
> $REDUX_ENV/new.1
#>> $REDUX_ENV/new.1

#this gives me the variants that haven't changed in a kit's data (according to the syntax of the variants.vcf file)
cat $REDUX_ENV/new.1 | grep "\.$" > $REDUX_ENV/old.0

#and this gives me the variants that have (again, according to the syntax)
cat $REDUX_ENV/new.1 | sed "/\.$/d" | sed -e "s/|/,/g"> $REDUX_ENV/new.2

#from the "new" variants, these are the ones that the hg38 list has names/definitions for
comm -13 $REDUX_ENV/hg38.0 $REDUX_ENV/new.2 > $REDUX_ENV/new.3

#to compare these "new variants" that hg38 knows about too ... let's just deal with the beginning locations (for the hg38 list)
cat $REDUX_ENV/hg38.0|cut -d',' -f1 | sort -u > $REDUX_ENV/hg38.1
#and same thing for the new list
cat $REDUX_ENV/new.3|cut -d',' -f1 | sort -u > $REDUX_ENV/new.4

#and then just isolate those where the new list is an even newer version of this variant
comm -12 hg38.1 new.4 |sed -e 's/$/,/g' > $REDUX_ENV/new.5
#hg38.1: 10147234,G,T
#new.4.: 10147234,G,A <-- truly new

#OK! now let's get the rest of the data (so we have the complete definition for the new variant)
grep -f $REDUX_ENV/new.5 $REDUX_ENV/hg38.0 > $REDUX_ENV/new.6

