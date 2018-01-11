#!/bin/python

# authors/licensing {{{

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,argparse,yaml,os,glob,shutil,re
import time,sqlite3,csv,zipfile,numpy as np
from collections import defaultdict
from lib import *
from db import *

# }}}

# trace

start_time = time.clock()
trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))

# env

try:
    REDUX_CONF = os.environ['REDUX_CONF']
except:
    print "Missing environment variable REDUX_CONF. Aborting."
    sys.exit()
try:
    REDUX_ENV = os.environ['REDUX_ENV']
except:
    print "Missing environment variable REDUX_ENV. Aborting."
    sys.exit()

config = yaml.load(open(REDUX_CONF))
#print ('\n'+sys.argv[0]+' version: '+config['VERSION']+'\n')

# arg parser (we can replace this -if useful- with getOpts or whatever)

parser = argparse.ArgumentParser()
parser.add_argument('-A', '--all', help='perform all possible steps', action='store_true')
parser.add_argument('-b', '--backup', help='do a backup', action='store_true')
parser.add_argument('-p', '--prep', help='prep file structure', action='store_true')
parser.add_argument('-d', '--data', help='SNP data processing', action='store_true')
args = parser.parse_args()

print ""

def go_backup():

    print "** performing backup."
    print "** (msg) CREATING BACKUP COPIES OF EXISTING FILES..."
    
    # autobackup dir
    refresh_dir('autobackup')
    for FILE_PATTERN in config['backup_files'].split():
        for FILE in glob.glob(FILE_PATTERN):
            shutil.copy(FILE,'autobackup')
    
    # autobackup2 dir {{{

    # Make further backup copies when running the script from scratch
    # This is useful when you want to make changes to the bad/inconsistent list, but still want to compare to the original previous run.
    # For example:
    # gawk 'NR==FNR {c[$5]++;next};c[$5]==0' tree.txt autobackup2/tree.txt
    # will tell you the changes to the tree structure that have resulted from the addition of new kits between "from-scratch" runs.
    
    #refresh_dir('autobackup2')
    #for FILE_PATTERN in config['backup_files'].split():
    #    for FILE in glob.glob(FILE_PATTERN):
    #        shutil.copy(FILE, 'autobackup2')
    
    # }}}
    
    if config['make_report']:
        #print "MAKING REPORT..."
        delete_file('report.csv')
    
    print "** + backup done."
    
def go_prep():

    print "** prepare file structure."

    # SKIPZIP check (beg)

    config['skip_zip'] = False
    if not config['skip_zip']:

        # Check ZIPDIR - contains existing zip files {{{

        if not os.path.exists(config['zip_dir']):
            print "Input zip folder does not appear to exist. Aborting.\n"
            sys.exit()

        # }}}
        # Check WORKING - the zip working exists && Empty it, otherwise make it {{{

        refresh_dir('working')

    # }}}
        # Check UNZIPDIR - contains zip output; refresh it {{{

        refresh_dir('unzips',not config['zip_update_only'])

        # }}}
        # Get the list of input files {{{

        FILES = glob.glob(REDUX_ENV+'/zips/bigy-*.zip')

        if len(FILES) == 0:
            print "No input files detected in zip folder. Aborting."
            print "Check name format: should be bigy-<NAME>-<NUMBER>.zip\n"
            sys.exit()
        else:
            print "input files detected: " + str(FILES)

        # }}}
        # Check whether unzip is installed {{{

        if not cmd_exists('unzip'):
            print "Unzip package not found. Aborting."
            sys.exit()

        # }}}
        # Check whether SNP list exists {{{

        csv = config['SNP_CSV']
        if not os.path.exists(csv):
            print "SNP names file does not exist. Try:"
            print "wget http://ybrowse.org/gbrowse2/gff/"+csv+" -O "+csv+"\n"
            sys.exit()

        # }}}
        # Check whether merge-ignore list exists {{{

        touch_file('merge-ignore.txt')

        # }}}
     
        # fix bash code (beg)

        # Unzip each zip in turn {{{

        print "Unzipping..."

        if config['zip_update_only']:
            #FILES=(`diff <(ls zip/bigy-*.zip | sed 's/zip\/bigy-//' | sed 's/.zip//') <(ls unzip/*.vcf | sed 's/unzip\///' | sed 's/.vcf//') | grep '<' | awk '{print "zip/bigy-"$2".zip"}'`)
            SET = [set(re.sub('zip/bigy-','',re.sub('.zip','',S)) for S in glob.glob('zips/bigy-*.zip'))]-set([re.sub('bigy-','',S) for S in glob.glob('unzips/*.vcf')])
            #print  ${#FILES[@]} "new files found"
            print "new files found: "+len(SET)
            print "new files detected: " + list(SET)

        #FILECOUNT=0

        #for ZIPFILE in ${FILES[@]}; do

        #    let FILECOUNT+=1
        #    PREFIX=`echo "$ZIPFILE" | gawk -F- '{print $2"-"$3}' | sed 's/.zip//'`
        #    #echo $FILECOUNT: $ZIPFILE : $PREFIX
        #    unzip -q $ZIPFILE -d working/
        #    if [ -s working/*.vcf ]; then mv working/*.vcf working/"$PREFIX".vcf; fi
        #    if [ -s working/*.bed ]; then mv working/*.bed working/"$PREFIX".bed; fi
        #    if [ -s working/*/variants.vcf ]; then mv working/*/variants.vcf working/"$PREFIX".vcf; fi
        #    if [ -s working/*/regions.bed ]; then mv working/*/regions.bed working/"$PREFIX".bed; fi
        #    if [ -s working/"$PREFIX".vcf ] && [ -s working/"$PREFIX".bed ]; then
        #        mv working/"$PREFIX".vcf unzip/;
	#	mv working/"$PREFIX".bed unzip/;
        #    else echo ""; echo "Warning: could not identify VCF and/or BED file for $PREFIX"
        #    fi

        #    rm -r working; mkdir working
        #    echo -n "."

        #done

        #echo ""

        # }}}

        # fix bash code (end)
    #fi

    # SKIPZIP check (end)

    # fix bash code (beg)

    # Skip some more if SKIPZIP set {{{

    #if config['skip_zip'] > 1:
    #    cp header.csv report.csv
    #    NFILES=`head -1 header.csv | gawk -v FS=, '{print NF-17}'`
    #    echo "... $NFILES results to be post-processed"

    #if config['skip_zip'] < 1:
    #    # Check number of BED = number of VCF files
    #    if [ `ls unzip/*.bed | wc -l` != `ls unzip/*.vcf | wc -l` ]; then
    #    echo "Number of BED files does not equal number of VCF files."
    #    echo "This is an unexpected error. Aborting."
    #    sys.exit()

    #T1=`date +%s.%N`
    #DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
    #echo "...complete after $DT seconds"

    # }}}
    # Generate statistics from BED & VCF files {{{

    #echo "Generating preliminary statistics..."

    #FILES=(`ls unzip/*.bed`)
    #echo "Total Kits:,,,,,"${#FILES[@]}',,,,,,,,,,,Kit' > header.csv
    #echo 'KEY:,,,,,,,,,,,,,,,,Date' >> header.csv
    #echo 'N+/N-,Number of +/- calls,,,,,,,,,,,,,,,Coverage' >> header.csv
    #echo '(?+),Call uncertain but presumed positive,(position forced),,,,,,,,,,,,,,...for age analysis' >> header.csv
    #echo 'cbl,Occurs on lower boundary of coverage,(often problematic),,,,,,,,,,,,,,Regions' >> header.csv
    #echo 'cbu,Occurs on upper boundary of coverage,(usually ok),,,,,,,,,,,,,,Variants' >> header.csv
    #echo 'cblu,Occurs as a 1-base-pair region,,,,,,,,,,,,,,,Passed' >> header.csv
    #echo '1stCol,First column which is positive,,,,,,,,,,,,,,,Simple SNPs' >> header.csv
    #echo 'Recur,Recurrencies in tree,(check: 1 or (R)),,,,,,,,,,,,,,SNPs under' "$TOPSNP" >> header.csv
    #echo '(s?),Questionable singleton,(not negative in some clademates),,,,,,,,,,,,,,Singleton SNPs' >> header.csv
    #echo '(s?!),Questionable singleton,(not negative in all clademates),,,,,,,,,,,,,,...for age analysis' >> header.csv
    #echo '(R),Allowed recurrency,,,,,,,,,,,,,,,Indels' >> header.csv
    #echo 'Blank,Securely called negative,,,,,,,,,,,,,,,Indels under' "$TOPSNP" >> header.csv
    #echo 'Full report at:,www.jb.man.ac.uk/~mcdonald/genetics/report.csv,,,,,,,,,,,,,,,Singleton Indels' >> header.csv
    #echo 'Non-shared SNPs' >> header.csv
    #echo 'GrCh37,Name(s),Ref,Alt,Type,N+,(?+),N-,nc,cbl+,cbl-,cbu+,cbu-,cblu+,cblu-,1stCol,Recur' >> header.csv
    #echo "Generating statistics for" ${#FILES[@]} "BED files..."
    #echo -n '[1/5] '
    #KITNAMES=`ls unzip/*.bed | sed 's/unzip\///g' | sed 's/.bed//g' | awk '1' ORS=,`
    #echo -n '[2/5] '
    #KITDATES=`ls -l --time-style +%Y-%m-%d unzip/*.bed | cut -d\  -f6 | awk '{print}' ORS=,`
    #echo -n '[3/5] '
    #STATS1=`gawk 'NR==FNR {a[NR]=$2;b[NR]=$3;n=NR} FNR==1 && NR!=1 {if (nfiles>0) print s,as,nrf;s=as=0;nfiles++} NR!=FNR {s+=$3-$2; for (i=1;i<=n;i++) if ($2<=b[i] && $3>=a[i]) {x=($3>b[i]?b[i]:$3)-($2>a[i]?$2:a[i]); if (x<0) x=0; as+=x}} {nrf=FNR} END {print s,as,FNR}' age.bed unzip/*.bed`
    #echo -n '[4/5] '
    #STATS2=`gawk '$1=="chrY" {n++} $1=="chrY" && $7=="PASS" {v++; if ($4!="." && $5!=".") {if (length($4)==1 && length($5)==1) {s++} else {i++} }} FNR==1 && NR!=1 {print n,v,s,0,0,0,i,0,0; n=v=s=i=0} END {print n,v,s,0,0,0,i,0,0}' unzip/*.vcf`

    #echo -n '[5/5] '
    #echo "$KITNAMES" | awk '{print substr($0,1,length($0)-1)}' > foo
    #echo "$KITDATES" | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS1" | awk '{print $1}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS1" | awk '{print $2}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS1" | awk '{print $3}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $1}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $2}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $3}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $4}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $5}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $6}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $7}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $8}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #echo "$STATS2" | awk '{print $9}' ORS=, | awk '{print substr($0,1,length($0)-1)}' >> foo
    #paste header.csv foo | sed 's/\t/,/' > fubar
    #mv fubar header.csv

    ## This does the same thing, but slower. From version 0.6.1
    ## for BEDFILE in ${FILES[@]}; do
    ##	VCFFILE=`echo "$BEDFILE" | sed 's/.bed/.vcf/'`
    ##	KITNAME=`echo "$BEDFILE" | gawk -F/ '{print $2}' | sed 's/.bed//'`
    ##	KITDATE=`ls -l --time-style +%Y-%m-%d "$BEDFILE" | cut -d\  -f6`
    ##	STATS=`gawk 'NR==FNR {a[NR]=$1;b[NR]=$2;n=NR} NR!=FNR {s+=$3-$2-; for (i=1;i<=n;i++) if ($2<=b[i] && $3>=a[i]) {x=($3>b[i]?b[i]:$3)-($2>a[i]?$2:a[i]); if (x<0) x=0; as+=x}} END {print s,as,FNR}' age.bed "$BEDFILE"`
    ##	STATS2=`gawk '$1=="chrY" {n++} $1=="chrY" && $7=="PASS" {v++; if ($4!="." && $5!=".") {if (length($4)==1 && length($5)==1) {s++} else {i++} }} END {print n,v,s,0,0,0,i,0,0}' "$VCFFILE"`
    ##	STATS="$KITNAME $KITDATE $STATS $STATS2"
    ##	gawk -v s="$STATS" 'NR==1 {split(s,stat," ")} {print $0","stat[NR]}' header.csv > foo
    ##	mv foo header.csv
    ##	echo -n "."
    ##done

    #echo ""
    #cp header.csv report.csv

    #T1=`date +%s.%N`
    #DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
    #echo "...complete after $DT seconds"

    # Close SKIPZIP if

    #fi

    # }}}
    # Skip some more if SKIPZIP set {{{

    #if [ "$SKIPZIP" == "0" ]; then

    # }}}
    # Identify list of variants {{{

    #echo "Identifying list of variants..."
    #gawk '$1=="chrY" && $7=="PASS" && $4!="." && $5!="." {print $2"\t"$4"\t"$5}' unzip/*.vcf | sed 's/,/;/g' > variant-list.txt

    #rm -f variant-list.txt
    #for BEDFILE in ${FILES[@]}; do
    #	VCFFILE=`echo "$BEDFILE" | sed 's/.bed/.vcf/'`
    #	gawk '$1=="chrY" && $7=="PASS" && $4!="." && $5!="." {print $2"\t"$4"\t"$5}' "$VCFFILE" | sed 's/,/;/g' >> variant-list.txt
    #	echo -n "."
    #done
    #echo ""

    # }}}
    # Add "missing" clades from file {{{

    # ! marks the implication so that is not counted when the SNP counts are made in the next section

    #gawk '$1=="^" {print $2"\t"$4"\t"$5"\t!"}' implications.txt >> variant-list.txt

    # }}}
    # Create a unique list of variants {{{

    #sort -nk1 variant-list.txt | uniq -c | sort -nk2 | gawk '{n="SNP"} $5=="!" {$1=0} length($3)>1 || length($4)>1 {n="Indel"} {print $2",,"$3","$4","n","$1",,,,,,,,,,,"}' > foo; mv foo variant-list.txt

    #T1=`date +%s.%N`
    #DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
    #print "...complete after $DT seconds"

    # }}}
    # Write out positive cases {{{

    #echo "Identifying positives and no calls..."

    # }}}
    # Include python script by Harald A. {{{

    #./positives-and-no-calls.py ${FILES[@]} > variant-match.txt

    #T1=`date +%s.%N`
    #DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
    #print "...complete after $DT seconds"

    # }}}

    # fix bash code (end)

    print "** + prep done."
    
def go_data():
    print "** process SNP data."
    print "** + SNP processing done."

if args.all:
    go_backup()
    go_prep()
    go_data()
else:
    if args.backup:
        go_backup()
    if args.prep:
        go_prep()
    if args.data:
        go_db()

print "** script complete.\n"

