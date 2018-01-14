#!/bin/bash

#Note: parse true new variants from an hg38 variants.vcf file

export TIMEFORMAT="time:%R"
echo;echo "step 0 - new.0 - iain's script, gets the variant data from new file"
echo "$ ./getVCFvariants.sh (iain's awk)> new.0"
time ./getVCFvariants.sh foo > new.0

cd $REDUX_ENV

echo;echo "step 1 - hg38.0 - get the default h38 data, doesn't need to be run every time"
echo "$ tail -n +2 $REDUX_ENV/snps_hg38.csv|cut -d'"' -f8,22,24 | sed -e 's/"/,/g' |sort -u > hg38.0"
time tail -n +2 $REDUX_ENV/snps_hg38.csv|cut -d'"' -f8,22,24 | sed -e 's/"/,/g' |sort -u> $REDUX_ENV/hg38.0

#echo "drop table variants_sh;" > $REDUX_ENV/go.tmp
#echo "create table variants_sh (allinfo char(20), code1 char(10), code2 char(10), variant int, tot tinyint, pass tinyint, fail tinyint);" \
#    >> $REDUX_ENV/go.tmp
echo;echo "step 2 - new.1 - boil down iain's data a bit more so I can use it for comparison" 
echo "$ [another awk command, mine] > new.1"
time cat $REDUX_ENV/new.0|awk 'BEGIN{p=0;f=0;t=0;a="beg";k="";c1="";c2="";}{
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

echo;echo "step 3 - old.0 - get variants that haven't changed (according to variants.vcf syntax)"
echo "$ cat new.1 | grep '\\.$' > old.0"
time cat $REDUX_ENV/new.1 | grep "\.$" > $REDUX_ENV/old.0

echo;echo "step 4 - new.2 - get variants that have changed (according to variants.vcf syntax)"
echo "$ cat new.1 | sed '/\\.$/d' | sed -e 's/|/,/g' > new.2"
time cat $REDUX_ENV/new.1 | sed "/\.$/d" | sed -e "s/|/,/g"> $REDUX_ENV/new.2

echo;echo "step 5 - new.3 - from the 'new' list, these are the ones hg38 has definitions for"
echo "$ comm -13 hg38.0 new.2 > new.3"
time comm -13 $REDUX_ENV/hg38.0 $REDUX_ENV/new.2 > $REDUX_ENV/new.3

echo;echo "step 6 - hg38.1 - let's boil hg38 reference to just the unique locations to compare"
echo "$ cat hg38.0|cut -d',' -f1 | sort -u > hg38.1"
time cat $REDUX_ENV/hg38.0|cut -d',' -f1 | sort -u > $REDUX_ENV/hg38.1
echo;echo "step 7 - new.4 - let's boil down new kit variant locations too to compare"
echo "$ cat new.3|cut -d',' -f1 | sort -u > new.4"
time cat $REDUX_ENV/new.3|cut -d',' -f1 | sort -u > $REDUX_ENV/new.4

echo;echo "step 8 - new.5 - now let's use comm to trap which one has new defs in the kit data"
echo "$ cat comm -12 hg38.1 new.4 |sed -e 's/$/,/g' > new.5"
time comm -12 hg38.1 new.4 |sed -e 's/$/,/g' > $REDUX_ENV/new.5
#hg38.1: 10147234,G,T
#new.4.: 10147234,G,A <-- truly new

echo;echo "step 9 - new.6 - and here's where we get the rest of the data for that new data"
echo "$ grep -f new.5 hg38.0 > new.6"
time grep -f $REDUX_ENV/new.5 $REDUX_ENV/hg38.0 > $REDUX_ENV/new.6

echo "step 10 - get count of unique new variants"
cntWInsDel=`cat $REDUX_ENV/new.6 | wc -l`
cntWOInsDel=`cat $REDUX_ENV/new.6 | sed '/ins/d' |wc -l`
echo "count w/ ins+del (bad entries): $cntWInsDel"
echo "count w/o: $cntWOInsDel"



