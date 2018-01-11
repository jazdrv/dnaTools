#!/bin/python

# {{{

@author: Iain McDonald
Contributors: Jef Treece, Harald Alvestrand

Purpose: Reduction and comparison script for Y-chromosome NGS test data
For free distribution under the terms of the GNU General Public License,
version 3 (29 June 2007)
https://www.gnu.org/licenses/gpl.html

# }}}

import sys,argparse,yaml,os,glob,shutil,re
import time,sqlite3,csv,shutil,zipfile,numpy as np
from collections import defaultdict

# routines - debug

def trace (level, msg):
    if level <= verbosity:
        print(msg)
    

# routines - file/dir - Zak

def refresh_dir(DIR,cleanFlag=False):
    DIR = REDUX_ENV+'/'+DIR
    if (os.path.isdir(DIR)):
        files = glob.glob(DIR+'/*')
        if cleanFlag:
            for f in files:
                os.remove(f)
    else:
        os.makedirs(DIR)
    
def delete_file(FILE):
    FILE = REDUX_ENV+'/'+FILE
    if os.path.exists(FILE):
        os.remove(FILE)
    
def touch_file(FILE):
    FILE = REDUX_ENV+'/'+FILE
    if not os.path.exists('merge-ignore.txt'):
        open('merge-ignore.txt','w').close()
    
def cmd_exists(CMD):
    return any(os.access(os.path.join(path, CMD), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

# routines - file/dir - Jef/Harald

def setup_dirs (unzip_dir):
    shutil.rmtree(unzip_dir, ignore_errors=True)
    os.makedirs(unzip_dir)
    
def extract_zips(unzip_dir,zip_dir):

    if not os.path.isdir(zip_dir):
        trace (0, '   Warn: no directory with zip files: %s' % zip_dir)
        return []

    FILES=os.listdir(zip_dir)

    # try to parse out at least the kit number by trying a series of regular expressions
    # adding regular expressions at the end of this list is safer than at the beginning
    # order is important - rules at top are matched first

    # constants used in filename regular expressions
    # groupings (?:xxx) are ignored
    ws = r'^[_]?'
    nam1 = r"[a-z]{0,20}|O\&#39;[a-z]{3,20}|O['][a-z]{3,20}"
    cname = r'([\w]{1,20})' #matches unicode chars; also matches digits though
    pnam = r'\('+nam1+r'\)'
    nam2 = r'(?:' +nam1 +'|' +pnam +r')'
    ndate = r'(?:(201[1-8][\d]{4}|201[1-8]-\d\d-\d\d|\d{4}201[1-8]))'
    sep = r'[\-\s\._]'
    seps = r'[\-\s\._]?'
    # sepp = r'[\-\s\._]+'
    sepp = r'_' # only use underscore as field separator
    sept = r'[\-\s\._]{3}'
    bigy = r'(?:big' +seps+ r'y(?:data)?|ydna)'
    rslt = r'(?:results|data|rawdata|vcfdata|raw data|csvexport|raw_data|raw|bigyrawdata)'
    name = r'((?:'+nam2+seps+'){1,3})'
    kit = r'(?:(?:kit|ftdna)?[ #]?)?([enhb1-9][0-9]{3,6})'
    rzip = r'zip(?:.zip)?'
    plac = r'([A-Z]{2})'
    name_re = [
        #0 e.g. bigy-Treece-N4826.zip
        (re.compile(ws+sep.join([bigy,name,kit,rzip]), re.I), 'name', 'kit'),
        #1 e.g. N4826_Treece_US_BigY_RawData_2018-01-03.zip
        (re.compile(ws +sepp.join([kit,name,plac,bigy,rslt,ndate])+'.zip', re.I), 'kit', 'name'),
        #2 e.g. 548872_Lindström_Germany_BigY_RawData_2018-01-01.zip
        (re.compile(ws +sepp.join([kit,cname,plac,bigy,rslt,ndate])+'.zip', re.I), 'kit', 'name'),
        ]


    trace (25, '   File names mapped, according to which regular expression:')
    # track counts - only for diagnostics
    cnt = defaultdict(int)
    # list of non-matching files
    nomatch=[]
    # all of the file names we could parse
    fname_dict = {}

    for line in FILES:
        fname = line.strip()
        if fname in rename_dict:
            # hand-edited filename mappings
            kkit, nname = rename_dict[fname]
            fname_dict[fname] = kkit, nname
            trace(25, '     {3:>2} {0:<50s}{1:<15s}{2:<10s}'.format(fname, nname, kkit, 'd'))
            cnt['d'] += 1
        else:
            if fname[-4:] not in ('.zip'):
                trace (15, '   Found foreigner hanging out in zip directory: {0}'.format(fname))
                continue
            d = {}
            for ii, (r,k1,k2) in enumerate(name_re):
                s = r.search(line)
                if s:
                    d[k1] = s.groups()[0]
                    if k2:
                        d[k2] = s.groups()[1]
                    else:
                        d['name'] = 'Unknown'
                    try:
                        trace (25, '     {3:>2} {0:<50s}{1:<15s}{2:<10s}'.format(fname,
                                                   d['name'], d['kit'], ii))
                        cnt[ii] += 1
                        fname_dict[fname] = d['kit'], d['name']
                    except:
                        trace (1, '   FAILURE on filename:', fname)
                    break
            else:
                nomatch.append(line)

    trace (20, '   Number of filenames not matched: {0}'.format(len(nomatch)))
    trace (22, '   Which expressions were matched:')
    for nn,cc in cnt.items():
        trace (22, '     {0:>2}: {1:>4}'.format(nn,cc))

    if len(nomatch) > 0:
        trace (10, '   Files that did not match:')
        for ll in nomatch:
            trace (10, '    %s' % ll.strip())

    # keep track of what needs to be cleaned up
    emptydirs = []

    for fname in fname_dict:
        kitnumber, kitname = fname_dict[fname]
        try:
            zf = zipfile.ZipFile(os.path.join(zip_dir, fname))
        except:
            trace (1, '   ERROR: file %s is not a zip' % fname)
        listfiles = zf.namelist()
        bedfile = vcffile = None
        for ff in listfiles:
            dirname, basename = os.path.split(ff)
            if basename == 'regions.bed':
                bedfile = ff
            elif basename == 'variants.vcf':
                vcffile = ff
            if dirname and (dirname not in emptydirs):
                emptydirs.append(dirname)
        if (not bedfile) or (not vcffile):
            trace(1, '   Warn: missing data in '+fname)
            continue
        if (bedfile == None) ^ (vcffile == None):
            trace(1, '   Warn: BED or VCF file is missing for %s' % fname)
            trace(1, '   This is an unexpected error. %s not processed.' % fname)
            continue
        zf.extractall(unzip_dir, [bedfile, vcffile])
        base = '%s-%s' % (kitname, kitnumber)
        try:
            fpath = os.path.join(unzip_dir, '%s')
            trace (40, "      "+fpath % base)
            os.rename(fpath % bedfile, (fpath % base)+'.bed')
            os.rename(fpath % vcffile, (fpath % base)+'.vcf')
        except:
            trace(1, '   Warn: could not identify VCF and/or BED file for '+base)

    # clean up any empty dirs unzip created
    if emptydirs:
        trace (30, '   Trying to remove droppings:')
        for dir in emptydirs:
            try:
                dp = os.path.join(unzip_dir, dir)
                os.removedirs(dp)
                trace (30, '     {0}'.format(dp))
            except FileNotFoundError:
                pass
            except:
                trace (30, '     W! could not remove {0}'.format(dp))
                pass

    # list of file names we unzipped
    files = os.listdir(unzip_dir)
    return files
def unpack(zip_dir, unzip_dir, verbosity):

    # messy problem - messy solution - kit names not consistent - Harald/Jef
    # collect run time statistics

    trace(10,'   Running the unpack-zip script...')
    setup_dirs(unzip_dir)
    fnames = extract_zips(unzip_dir, zip_dir)
    trace (10, '   Number of files: {0}'.format(len(fnames)))
    trace (40, '   Files unpacked:')
    for ff in fnames:
        trace (40, ff)

# routines - arghandler - Zak

def go_all():
    print "go all"
    # go_backup()
    # go_prep()
    # go_data()
    
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

# SNP extraction routines based on original - Harald 
# extracts the SNP calls from the VCF files and
# determines the coverage of SNPs in the BED files of BigY tests.

def analyzeVcf(file):

    #Returns a dict of position -> mutation mappings

    with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
        trace (30, "   Extracting VCF: %s" % vcffile)
        result = {}
        for line in vcffile:
            fields = line.split()
            if (fields[0] == 'chrY' and fields[6] == 'PASS' and fields[3] != '.' and fields[4] != '.'):
                # fix by Jef Treece for fields containing commas:
                result[int(fields[1])] = fields[1] + '.' + fields[3].replace(',', ';') + '.' + fields[4].replace(',', ';')
                # result[int(fields[1])] = fields[1] + '.' + fields[3] + '.' + fields[4]
        return result
    
def analyzeBed(file):

    #Returns an array of path segments.

    with open(os.path.splitext(file)[0] + '.bed') as bedfile:
        trace (30, "   Extracting BED: %s" % bedfile)
        result = []
        for line in bedfile:
        fields = line.split()
        if (fields[0] == 'chrY'):
            result.append((int(fields[1]), int(fields[2])))
    return result
    
def makeCall(pos, index_container, bed_calls):

    #Figure out whether this position is on a segment boundary.
    #Between segments = 'nc'; top of segment = 'cbu'; bottom of segment = 'cbl'.
    #Only call in a single-position segment = 'cblu'.
    #index_container contains first segment to be looked at.
    #This function must only be called for increasing values of pos, and with
    #sorted bed_calls.

    call = ';nc'
    for bed_index in range(index_container[0], len(bed_calls)):
        pos_pair = bed_calls[bed_index]
        index_container[0] = bed_index
        if pos_pair[1] >= pos:
            # Position is before or within this segment.
            if pos_pair[0] <= pos:
                # Position is within this segment.
                if pos_pair[0] == pos_pair[1] and pos_pair[0] == pos:
                    call = ';cblu'
                elif pos_pair[0] == pos:
                    call = ';cbl'
            elif pos_pair[1] == pos:
                call = ';cbu'
            else:
                call = ''
        else:
            # Position is before this segment.
            call = ';nc'
        return call
        # If position is after segment, continue.
    return ';nc' # After end of last segment.
    
def extract(unzip_dir,files,variants):

  d = []
  s = []

  curpath = os.path.abspath(os.curdir)
  with open(os.path.join(curpath, 'variant-list.txt')) as line_headings:
    for line in line_headings:
      d.append(line.rstrip())
      x = line.split(',')
      s.append(int(x[0]))  # s holds the genome position for each line

  for file in files:
    vcf_calls = analyzeVcf(unzip_dir + file)
    bed_calls = analyzeBed(unzip_dir + file)
    bed_index = [0]
    for lineno in range(len(d)):
      d[lineno] += ','
      if s[lineno] in vcf_calls:
        d[lineno] += vcf_calls[s[lineno]]
      d[lineno] += makeCall(s[lineno], bed_index, bed_calls)

  for line in d:
    print (line)
def file_len(fname):

    #File length, thanks to StackOverflow
    #https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python

    i=-1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# routines - Iain 

def readVcf(file):

    #Returns a dict of position -> mutation mappings
    #Modified from Harald's analyzeVCF, this version returns every mutation with
    #its derived value, regardless of whether it was ancestral or not

    with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
        trace (30, "   Extracting VCF: %s" % vcffile)
        result = {}
        for line in vcffile:
            fields = line.split()
            if (fields[0] == 'chrY' and int(fields[1]) > 0 and fields[3] != '.' and fields[4] != '.'):
                result[fields[1]] = [int(fields[1]), str(fields[3]), str(fields[4])]
        return result

# Main

def zak():

    try:
        REDUX_CONF = os.environ['REDUX_CONF']
    except:
        print "Missing environment variable REDUX_CONF. Aborting."
        syx.exit()
    try:
        REDUX_ENV = os.environ['REDUX_ENV']
    except:
        print "Missing environment variable REDUX_ENV. Aborting."
        syx.exit()

    config = yaml.load(open(REDUX_CONF))

    parser = argparse.ArgumentParser()
    parser.add_argument('-A', '--all', help='perform all possible steps', action='store_true')
    parser.add_argument('-b', '--backup', help='do a backup', action='store_true')
    parser.add_argument('-p', '--prep', help='prep file structure', action='store_true')
    parser.add_argument('-d', '--data', help='SNP data processing', action='store_true')
    args = parser.parse_args()

    #print ('\n'+sys.argv[0]+' version: '+config['VERSION']+'\n')
    print ""

    if args.all:
        go_all()
    else:
        if args.backup:
            go_backup()
        if args.prep:
            go_prep()
        #if args.data:
            # go_data()

    print "** script complete.\n"
    
def main():

    # Set start time for efficiency testing {{{

    start_time = time.clock()
    trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))
    
    #}}}
    # Init variant database {{{

    trace (1, "Initialising database...")
    db = sqlite3.connect('variant.db')
    dc = db.cursor()

    #}}}
    # Drop table{{{

    if (droptables > 0):
        dc.execute('''DROP TABLE IF EXISTS variants''')
        dc.execute('''DROP TABLE IF EXISTS hg19''')
        dc.execute('''DROP TABLE IF EXISTS hg38''')
        dc.execute('''DROP TABLE IF EXISTS kits''')
        dc.execute('''DROP TABLE IF EXISTS people''')
        dc.execute('''DROP TABLE IF EXISTS strs''')
        dc.execute('''DROP TABLE IF EXISTS calls''')
        dc.execute('''DROP TABLE IF EXISTS tree''')

    #}}}
    # Create table of variants if it doesn't already exist{{{

    dc.execute('''CREATE TABLE IF NOT EXISTS variants
                  (id INTEGER PRIMARY KEY,
                  grch37 INTEGER,
                  length SMALLINT,
                  ref TEXT,
                  alt TEXT,
                  inage BOOLEAN)''')

    #}}}
    # Create reference tables of variants{{{

    dc.execute('''CREATE TABLE IF NOT EXISTS hg19
                  (id INTEGER PRIMARY KEY,
                  snp BOOLEAN,
                  grch37 INTEGER,
                  grch37end INTEGER,
                  name CHARACTER(32),
                  anc CHARACTER(64),
                  der CHARACTER(64))''')
    dc.execute('''CREATE TABLE IF NOT EXISTS hg38
                  (id INTEGER PRIMARY KEY,
                  snp BOOLEAN,
                  grch38 INTEGER,
                  grch38end INTEGER,
                  name CHARACTER(32),
                  anc CHARACTER(64),
                  der CHARACTER(64))''')

    #}}}
    # Create reference tables of kits{{{

    dc.execute('''CREATE TABLE IF NOT EXISTS kits
                  (id INTEGER PRIMARY KEY,
                  vcffile CHARACTER(256),
                  bedfile CHARACTER(256),
                  company CHARACTER(16),
                  test CHARACTER(16),
                  kitnum CHARACTER(16),
                  testdate INTEGER,
                  uploaddate INTEGER,
                  personid INTEGER)''')

    # }}}
    # Create reference tables of people{{{

    dc.execute('''CREATE TABLE IF NOT EXISTS people
                  (person INTEGER PRIMARY KEY,
                  ngstest BOOLEAN,
                  ftdnaid CHARAACTER(16) REFERENCES kits(kitnum),
                  fgcid CHARAACTER(16) REFERENCES kits(kitnum),
                  yseqid CHARAACTER(16) REFERENCES kits(kitnum),
                  otherid CHARAACTER(16) REFERENCES kits(kitnum),
                  mdkaname CHARACTER(16),
                  mdkacountry CHARACTER(2),
                  mdkalat REAL,
                  mdkalong REAL,
                  mdkalocweight REAL,
                  mdkadesc CHARACTER(255),
                  mdkayear SMALLINT,
                  mdkaunc SMALLINT,
                  mdkaafter SMALLINT,
                  mdkabefore SMALLINT,
                  mrkhid INTEGER REFERENCES tree(id),
                  coverage INTEGER,
                  agecoverage INTEGER,
                  totcoverage INTEGER,
                  nsnps INTEGER,
                  snpsintree SMALLINT,
                  singletons SMALLINT,
                  qsingletons SMALLINT,
                  agesingletons SMALLINT)''')

    #}}}
    # Create reference tables of STRs{{{

    dc.execute('''CREATE TABLE IF NOT EXISTS strs
                  (person INTEGER PRIMARY KEY,
                  DYS393 TINYINT,
                  DYS390 TINYINT,
                  DYS19 TINYINT,
                  DYS391 TINYINT,
                  DYS385a TINYINT,
                  DYS385b TINYINT,
                  DYS426 TINYINT,
                  DYS388 TINYINT,
                  DYS439 TINYINT,
                  DYS389i TINYINT,
                  DYS392 TINYINT,
                  DYS389ii TINYINT,
                  DYS458 TINYINT,
                  DYS459a TINYINT,
                  DYS459b TINYINT,
                  DYS455 TINYINT,
                  DYS454 TINYINT,
                  DYS447 TINYINT,
                  DYS437 TINYINT,
                  DYS448 TINYINT,
                  DYS449 TINYINT,
                  DYS464a TINYINT,
                  DYS464b TINYINT,
                  DYS464c TINYINT,
                  DYS464d TINYINT,
                  DYS460 TINYINT,
                  YH4 TINYINT,
                  YCAIIa TINYINT,
                  YCAIIb TINYINT,
                  DYS456 TINYINT,
                  DYS607 TINYINT,
                  DYS576 TINYINT,
                  DYS570 TINYINT,
                  CDYa TINYINT,
                  CDYb TINYINT,
                  DYS442 TINYINT,
                  DYS438 TINYINT,
                  DYS531 TINYINT,
                  DYS578 TINYINT,
                  DYF395S1a TINYINT,
                  DYF395S1b TINYINT,
                  DYS590 TINYINT,
                  DYS537 TINYINT,
                  DYS641 TINYINT,
                  DYS472 TINYINT,
                  DYF406S1 TINYINT,
                  DYS511 TINYINT,
                  DYS425 TINYINT,
                  DYS413a TINYINT,
                  DYS413b TINYINT,
                  DYS557 TINYINT,
                  DYS594 TINYINT,
                  DYS436 TINYINT,
                  DYS490 TINYINT,
                  DYS534 TINYINT,
                  DYS450 TINYINT,
                  DYS444 TINYINT,
                  DYS481 TINYINT,
                  DYS520 TINYINT,
                  DYS446 TINYINT,
                  DYS617 TINYINT,
                  DYS568 TINYINT,
                  DYS487 TINYINT,
                  DYS572 TINYINT,
                  DYS640 TINYINT,
                  DYS492 TINYINT,
                  DYS565 TINYINT,
                  DYS710 TINYINT,
                  DYS485 TINYINT,
                  DYS632 TINYINT,
                  DYS495 TINYINT,
                  DYS540 TINYINT,
                  DYS714 TINYINT,
                  DYS716 TINYINT,
                  DYS717 TINYINT,
                  DYS505 TINYINT,
                  DYS556 TINYINT,
                  DYS549 TINYINT,
                  DYS589 TINYINT,
                  DYS522 TINYINT,
                  DYS494 TINYINT,
                  DYS533 TINYINT,
                  DYS636 TINYINT,
                  DYS575 TINYINT,
                  DYS638 TINYINT,
                  DYS462 TINYINT,
                  DYS452 TINYINT,
                  DYS445 TINYINT,
                  YA10 TINYINT,
                  DYS463 TINYINT,
                  DYS441 TINYINT,
                  Y1B07 TINYINT,
                  DYS525 TINYINT,
                  DYS712 TINYINT,
                  DYS593 TINYINT,
                  DYS650 TINYINT,
                  DYS532 TINYINT,
                  DYS715 TINYINT,
                  DYS504 TINYINT,
                  DYS513 TINYINT,
                  DYS561 TINYINT,
                  DYS552 TINYINT,
                  DYS726 TINYINT,
                  DYS635 TINYINT,
                  DYS587 TINYINT,
                  DYS643 TINYINT,
                  DYS497 TINYINT,
                  DYS510 TINYINT,
                  DYS434 TINYINT,
                  DYS461 TINYINT,
                  DYS435 TINYINT)''')

    #}}}
    # Create reference matrix of calls{{{

    dc.execute('''CREATE TABLE IF NOT EXISTS calls
                  (id LONGINT PRIMARY KEY,
                  person INTEGER REFERENCES people(person),
                  variant INTEGER REFERENCES variants(id),
                  assigned BOOLEAN,
                  positive BOOLEAN,
                  callable BOOLEAN,
                  bq TINYINT,
                  mq TINYINT,
                  nder TINYINT,
                  nanc TINYINT
                  )''')

    #}}}
    # Create tree table{{{

    dc.execute('''CREATE TABLE IF NOT EXISTS tree
                  (id INTEGER PRIMARY KEY,
                  parendid INTEGER,
                  clade CHARACTER(16),
                  variants BLOB,
                  qualities BLOB,
                  ageflags BLOB,
                  children BLOB,
                  ancestralstrs BLOB,
                  originlat REAL,
                  originlong REAL,
                  coverage INTEGER,
                  agecoverage INTEGER,
                  poznikcoverage INTEGER,
                  combbedcoverage INTEGER,
                  olderthan SMALLINT,
                  olderthanunc SMALLINT,
                  olderthankit INTEGER,
                  youngerthan SMALLINT,
                  yongerthanunc SMALLINT,
                  youngerthankit1 INTEGER,
                  youngerthankit2 INTEGER,
                  snpage SMALLINT,
                  snpagelo SMALLINT,
                  snpagehi SMALLINT,
                  snpagepdf BLOB,
                  snpparentpdf BLOB,
                  strage SMALLINT,
                  stragelo SMALLINT,
                  stragehi SMALLINT,
                  stragepdf BLOB,
                  strparentpdf BLOB,
                  combage SMALLINT,
                  combagelo SMALLINT,
                  combagehi SMALLINT,
                  combagepdf BLOB)''')


    trace (1, "Processing Build 38 BigY files...")
    
    #}}}
    # Unpack ZIP files{{{

    if (skipto <= 1):
        trace (2, "Unpacking ZIP files...")
        unpack(zip_dir,unzip_dir,verbosity)
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)

        trace (5, "Associating unzipped files with kits...")
        
    #}}}
    # Associate kits with people{{{

    if (skipto <= 10):
        trace (2, "Associating kits with people...")
        
    #}}}
    # Generate dictionary of variant positions {{{

    #   The VCF files are parsed twice. The first pass identifies the list
    #   of variants to be queried. The second pass reads the calls for those
    #   variants. This means we only need to treat positions with a variant,
    #   rather than every position in a chromosome.
    #   Using a dictionary here means only one copy of each variant is
    #   created.

    if (skipto <= 11):

        trace (2, "Generating database of all variants...")
        vcffiles = [f for f in os.listdir(unzip_dir) if f.endswith('.vcf')]
        trace (10, "   %i files detected" % len(vcffiles))
        
        variant_dict = {}
        for file in vcffiles:
            vcf_calls = readVcf(unzip_dir + file)
            variant_dict.update(vcf_calls)
        
        trace (10, "   %i variants found" % len(variant_dict))
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)

        # Dump variant dictionary into sorted array {{{

        trace (20, "   Dumping variants into array...")
        variant_array = np.array(list(variant_dict.values()))

        #        variant_array = np.array([],dtype={'names': ('start', 'anc', 'der'),
        #                                          'formats': ('i4', 'S20', 'S20')})

        trace (30, "      Check variant [0] is %s" % variant_array[0])
        trace (30, "      Check variant [0] position is %s" % variant_array[0][1])
        trace (30, "      Check variant [%s] is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1]))
        trace (30, "      Check variant [%s] position is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1][1]))

        trace (20, "   Inserting data into variant array database...")
        dc.executemany('''INSERT INTO variants(id,ref,alt) VALUES (?,?,?)''', variant_array)

        #        Test data has entered database correctly
        #        dc.execute('SELECT * FROM variants LIMIT 5')
        #        print (dc.fetchone())

        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)

        }}}

    # }}}

    # Reading calls for variants {{{

    if (skipto <= 12):
        trace (2, "Generating database of calls...")
        vcffiles = [f for f in os.listdir(unzip_dir) if f.endswith('.vcf')]
        trace (10, "   %i files detected" % len(vcffiles))

        dc.execute('''INSERT INTO calls(variant,person)
                      SELECT id, person
                      FROM variants CROSS JOIN people''')
        # !!!


    #        Test data has entered database correctly
    #        dc.execute('SELECT * FROM calls LIMIT 5')
    #        print (dc.fetchone())

    # }}}
    # Name variants and derive ancestral values {{{

    # Some variants are positive in the reference sequence, so we need to
    # look up their ancestral values. We'll get the SNP names while we're
    # at it.

    if (skipto <= 13):
        trace (2, "Getting names of variants...")

        # Read in SNPs from reference lists
        trace (10, "   Importing SNP reference lists...")

        snp_reference = csv.reader(open(b37_snp_file))
        dc.executemany("INSERT INTO hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)",
                       ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
    
        snp_reference = csv.reader(open(b38_snp_file))
        dc.executemany("INSERT INTO hg38(grch38,grch38end,name,anc,der) VALUES (?,?,?,?,?)",
                       ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))

        # Test data has entered database correctly
        # dc.execute('SELECT * FROM hg38 LIMIT 5')
        # print (dc.fetchone())

        # Read in SNPs from reference lists
        # Probably doesn't need done at this point
        # trace (10, "   Joining reference lists to variant database...")


        # dc.execute('''SELECT hg38.grch38, hg38.name
        # FROM hg38
        # INNER JOIN hg19 on hg19.name = hg38.name''')

        # dc.execute('''SELECT variants.id, hg38.name
        # FROM variants
        # LEFT OUTER JOIN hg38 on hg38.grch38 = variants.id''')

        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
            
    # }}}
    # Print final message and exit {{{{

    t = float((time.clock() - start_time))
    trace (1, "Execution finished in: %.3f seconds" % t)

    # }}}

    # sys.exit(0)

main()

