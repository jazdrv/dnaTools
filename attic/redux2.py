#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 17:07:57 2017

@author: Iain McDonald

Contributors: Jef Treece, Harald Alvestrand, Zak Jones

Purpose: Reduction and comparison script for Y-chromosome NGS test data
For free distribution under the terms of the GNU General Public License,
version 3 (29 June 2007)
https://www.gnu.org/licenses/gpl.html

#{{{
*** Search below for USER DEFINED VARIABLES ***

Works on Anaconda Python
Windows Bash needs Numpy:
https://scipy.org/install.html
"""

# The following packages are required by this programme
import time #, psutil
#import numpy.lib.recfunctions as rfn
import csv, os
#import sys
import shutil
import zipfile
import re
import numpy as np
#from time import time
#from getopt import getopt
from collections import defaultdict
#from subprocess import call
from db import DB

# Current version
version="2.0.0.20180108.ALPHA"

# =============================================================================
# USER DEFINED VARIABLES

# This label will define what your top-level SNP is called out of the
# (often large number of) possible shared SNPs.
haplogroup="R"
topsnp="U106"

# This allows you to start the programme at a particular point:
    # skipto = 0 -> run whole programme
    # skipto = 10 -> work only with already-unzipped files
    # skipto = 20 -> do not re-extract the VCF calls
    # skipto = 30 -> do not re-make the tree
skipto=11

# This sets which steps you wnat to carry out [1] and which you don't [0]
makereport=1
makeages=1
makeshortreport=1
makehtmlreport=1

# This sets the mutation rate used to calculate the ages
# and it's 95% confidence interval
rate=0.8124 # x 10^-10 per base pair per year
rate0=0.7534 # x 10^-10 per base pair per year
rate1=0.8722 # x 10^-10 per base pair per year

# This sets what defines the present day
zeroage=1950.33   # AD
errzeroage=15.02  # years (standard deviation)

# Set verbosity of messages
verbosity = 30

# Sets directories and other input files
zip_dir = "zip/"
unzip_dir = "unzip/"

# The list of SNPs used as a reference
# Get this using:
#call(["source","webupdate.scr"])
# which runs:
#   wget http://ybrowse.org/gbrowse2/gff/snps_hg19.csv -O snps_hg19.csv
#   awk '{print $4,$5,$9,$11,$12}' FS='","' snps_hg19.csv | sort -nk1 | \
#   awk 'NR>1 {if (a!=$1) {print a,b,c,d,e,f; c=""}} {x=$0;a=$1;b=$2;c=length(c)>0?c"/"$3:$3;d=$4;e=$5;f=$6}' OFS='\t' > snps_hg19.tsv
# to convert it to TSV format with one line per position. Unfortunately,
# there doesn't seem to be a simple Pythonic way of reading in the original CSV
# file due to the quotes.
b37_snp_file = "snps_hg19.csv"
b38_snp_file = "snps_hg38.csv"

# This defines whether tables get made again
droptables = 1

# =============================================================================
# USER DEFINED FILE MAPPING
names = """
FTDNA345238Newell.zip, 345238, Newell
155941_BigY_RawData_20140911-1.zip, 155941, Unknown
Lee 237414 BigY Raw Data.zip, 237414, Lee
U106_515653_Hogenmiller_BigY_RawData_2016_11_20.zip, 515653, Hogenmiller
bigy-Bettinger57020.zip, 57020, Bettinger
"""
#-----end hand-editing
#}}}

rename_dict = {}

for row in csv.reader(names.splitlines()):
    if row and row[0]:
        rename_dict[row[0].strip()] = (row[1].strip(), row[2].strip())

# =============================================================================
# Here be dragons...
# =============================================================================
# SUBROUTINE DEFINIITIONS

# -----------------------------------------------------------------------------

"""
  Unpacking routines based on original by Jef Treece, 14 Mar 2017

  Purpose:
    This set of subroutines unpacks the zip files containing the
    VCF and BED files from BigY tests.
    
  Description:
    The input list of files from users are inhomogenous. In order to
    preserve a list of files linkable back to the original, a database is
    created in the programme to map file names to users' surnames and kit
    numbers.

  Original copyright:
    For free distribution under the terms of the
    GNU General Public License, version 3 (29 June 2007)
    https://www.gnu.org/licenses/gpl.html
"""

# -----------------------------------------------------------------------------
# trace (print) output if it exceeds the noise level
def trace (level, msg):
    if level <= verbosity:
        print(msg)

# -----------------------------------------------------------------------------
# create new directory for output
def setup_dirs (unzip_dir):
    shutil.rmtree(unzip_dir, ignore_errors=True)
    os.makedirs(unzip_dir)

# -----------------------------------------------------------------------------
# extract zip files
# messy problem - messy solution - kit names are not consistent
def extract_zips(unzip_dir, zip_dir):

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

# -----------------------------------------------------------------------------
def unpack(zip_dir, unzip_dir, verbosity):

    # collect run time statistics
    trace(10,'   Running the unpack-zip script...')
    setup_dirs(unzip_dir)
    fnames = extract_zips(unzip_dir, zip_dir)
    trace (10, '   Number of files: {0}'.format(len(fnames)))
    trace (40, '   Files unpacked:')
    for ff in fnames:
        trace (40, ff)

# -----------------------------------------------------------------------------

"""
  SNP extraction routines based on original by Harald Alvestrand, 2016

  Purpose:
    This set of subroutines extracts the SNP calls from the VCF files and
    determines the coverage of SNPs in the BED files of BigY tests.
    
  Implied copyright:
    For free distribution under the terms of the
    GNU General Public License, version 3 (29 June 2007)
    https://www.gnu.org/licenses/gpl.html
"""

# -----------------------------------------------------------------------------

def analyzeVcf(file):
  """Returns a dict of position -> mutation mappings"""
  with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
    trace (30, "   Extracting VCF: %s" % vcffile)
    result = {}
    for line in vcffile:
      fields = line.split()
      if (fields[0] == 'chrY' and fields[6] == 'PASS' and fields[3] != '.'
          and fields[4] != '.'):
# Fix by Jef Treece for fields containing commas:
        result[int(fields[1])] = fields[1] + '.' + fields[3].replace(',', ';') + '.' + fields[4].replace(',', ';')
#        result[int(fields[1])] = fields[1] + '.' + fields[3] + '.' + fields[4]
  return result


# -----------------------------------------------------------------------------

def analyzeBed(file):
  """Returns an array of path segments."""
  with open(os.path.splitext(file)[0] + '.bed') as bedfile:
    trace (30, "   Extracting BED: %s" % bedfile)
    result = []
    for line in bedfile:
      fields = line.split()
      if (fields[0] == 'chrY'):
        result.append((int(fields[1]), int(fields[2])))
  return result


# -----------------------------------------------------------------------------

def makeCall(pos, index_container, bed_calls):
  """Figure out whether this position is on a segment boundary.

  Between segments = 'nc'; top of segment = 'cbu'; bottom of segment = 'cbl'.
  Only call in a single-position segment = 'cblu'.
  index_container contains first segment to be looked at.
  This function must only be called for increasing values of pos, and with
  sorted bed_calls."""

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

# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# File length, thanks to StackOverflow
#https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python

def file_len(fname):
    i=-1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# -----------------------------------------------------------------------------

"""
  Other routines by Iain McDonald
"""

# -----------------------------------------------------------------------------

def readVcf(file):
  """
  Returns a dict of position -> mutation mappings
  Modified from Harald's analyzeVCF, this version returns every mutation with
  its derived value, regardless of whether it was ancestral or not
  """
  with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
    trace (30, "   Extracting VCF: %s" % vcffile)
    result = {}
    for line in vcffile:
      fields = line.split()
      if (fields[0] == 'chrY' and int(fields[1]) > 0 and fields[3] != '.' and fields[4] != '.'):
           result[fields[1]] = [int(fields[1]), str(fields[3]), str(fields[4])]
  return result

# =============================================================================
# MAIN PROGRAMME

def main():
    # Set start time for efficiency testing
    start_time = time.clock()
    trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))
    
    # -------------------------------------------------------------------------
    # Initialise variant database
    trace (1, "Initialising database...")
    db = DB()
    dc = db.dc

    # Drop table
    # Create reference tables of variants
    # Create reference tables of kits
    # Create reference tables of people
    # Create reference tables of STRs
    # Create reference matrix of calls
    # Create tree table
    # now handled by create_schema

    # =========================================================================
    trace (1, "Processing Build 38 BigY files...")
    
    # -------------------------------------------------------------------------
    # Unpack ZIP files
    if (skipto <= 1):
        trace (2, "Unpacking ZIP files...")
        unpack(zip_dir,unzip_dir,verbosity)
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)

        trace (5, "Associating unzipped files with kits...")
        
    # -------------------------------------------------------------------------
    # Associate kits with people
    if (skipto <= 10):
        trace (2, "Associating kits with people...")
        
    # -------------------------------------------------------------------------
    # Generate dictionary of variant positions
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

        # Dump variant dictionary into sorted array
        trace (20, "   Dumping variants into array...")
        variant_array = np.array(list(variant_dict.values()))
#        variant_array = np.array([],dtype={'names': ('start', 'anc', 'der'),
#                                          'formats': ('i4', 'S20', 'S20')})

        trace (30, "      Check variant [0] is %s" % variant_array[0])
        trace (30, "      Check variant [0] position is %s" % variant_array[0][1])
        trace (30, "      Check variant [%s] is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1]))
        trace (30, "      Check variant [%s] position is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1][1]))

        trace (20, "   Inserting data into variant array database...")
        dc.executemany('''INSERT INTO variants(id,anc,der) VALUES (?,?,?)''', variant_array)

#        Test data has entered database correctly
#        dc.execute('SELECT * FROM variants LIMIT 5')
#        print (dc.fetchone())

        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)

    # -------------------------------------------------------------------------
    # Reading calls for variants
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

    # -------------------------------------------------------------------------
    # Name variants and derive ancestral values
    #   Some variants are positive in the reference sequence, so we need to
    #   look up their ancestral values. We'll get the SNP names while we're
    #   at it.
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

#        Test data has entered database correctly
#        dc.execute('SELECT * FROM hg38 LIMIT 5')
#        print (dc.fetchone())

        # Read in SNPs from reference lists
#           Probably doesn't need done at this point
#        trace (10, "   Joining reference lists to variant database...")


#        dc.execute('''SELECT hg38.grch38, hg38.name
#                      FROM hg38
#                      INNER JOIN hg19 on hg19.name = hg38.name''')

#        dc.execute('''SELECT variants.id, hg38.name
#                      FROM variants
#                      LEFT OUTER JOIN hg38 on hg38.grch38 = variants.id''')

        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
            
    # Print final message and exit
    t = float((time.clock() - start_time))
    trace (1, "Execution finished in: %.3f seconds" % t)
#    sys.exit(0)



# =============================================================================
# RUN MAIN PROGRAMME
main()
