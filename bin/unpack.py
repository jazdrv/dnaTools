#!/usr/bin/env python3
# coding: utf-8
# The preference is python3; not well tested w/ python2
"""
  Purpose:
    * unpacks zip files containing vcf and bed files from BigY tests
    * corrects for inconsistent naming of data files
    * detects mistakes such as duplicate files
    * maps-out (ignores) known-bad files
    * provides a method to update/fix certain files
    * corrects (some) improperly packed .zip files
    * ability to utilize a partial subset of the .zip files

  Usage:
    unpack-zip-files [-z <zipdir>]+ -u <unzipdir> [-v*] [--dryrun]
    example: unpack-zip-files -z zip1 -z zip2 -u unzip -vvv > log.out

    Other optional command-line args:
    -k : keep already-unpacked zip files instead of removing them
    -r : copy and rename the zip files instead of extracting them
    -m <csvfile> : provide a filename mapping replacing default filemap.csv
    -s <bins>:<binnum> : only use a subset of the zip files; e.g. -s 2:1
       means split zips into two bins and use the first bin.
       binnum can also be a comma separated list; e.g. -s 10:1,3,5
       Special case: -s n:0 means use a random set that is 1/nth the total

    Process for handling new zip file(s):
    1. Put your zip files in zip.new directory
    2. Run unpack-zip-files.py -z zip.new -u unzip.new -vv -d > log.out
       (this is a dry run and won't save anything)
    3. Inspect log.out for warnings or errors. Usually obvious.
    4. Look at the list of files that didn't get mapped.
    5. If there are a large number, you could update regular expressions in
       the main program.
    6. More likely, you'll just add new lines to the csv spreadsheet file
    7. Repeat this until all files are mapped
    8. Run without the -d flag on final run

  Setup:
    Requires mapping spreadsheet (default name filemap.csv)
    General instructions are below.

  Copyright:
    For free distribution under the terms of the
    GNU General Public License, version 3 (29 June 2007)
    https://www.gnu.org/licenses/gpl.html

  Jef Treece
  v1.15, 3 Jan 2018


Problem: .zip files are saved with various naming standards across projects
    that utilize Big-Y data files.  A researcher needs to extract these files
    in a standard way.

Complications: Contributors used various naming standards. People have picked
    up these files, re-named them to their own personal preference, and
    perhaps even unzipped/re-zipped them, changing the original file
    signature. Furthermore, some files have been uploaded with errors,
    replaced with newer versions, and various other issues. The projects are
    distributed, and there is no central file repository, so globally
    re-naming files at the source is currently out of the question. Original
    owners of the data files sometimes want to make corrections to surname or
    kit number or even bad data that they would hope get pushed out to the
    various people who might be publishing results based on the files.

Implications: without standardizing, no researcher can build on another
    researcher's work without entirely depending on that person to supply the
    data. This is not scalable, and it's bottle-necked on correcting
    mistakes. At best, each researcher must re-do the work of re-naming. At
    worst, no results are reproducible and work is lost if someone disappears
    or is unresponsive. Without codifying this, everyone will have trouble
    correcting mistakes and accepting updates from the user base.

Strategy: This program aims to address the above issues. Firstly, the regular
    expressions in the main program do the heavy lifting by parsing out the
    surname and kit number from the input data file. If the file is named in
    any reasonable fashion, this is possible. If parsing is not possible,
    here, we handle exceptions and correct mistakes. By making this an open
    source script, people can supply corrections, and only the script need be
    updated and made available to anyone who wants to contribute. It is
    practical for multiple people to contribute. Additionally, this script
    strives to work with the file in its original form and considers it
    read-only. ORIGINAL DATA FILES should not be changed or re-named unless
    absolutely necessary (or as part of a planned and coordinated
    cross-project effort)! There are quite a few files in here that have been
    renamed by someone already, and that creates duplicate lines to maintain
    here. Retain the different variants, all map to the same output.

Directions for adjusting infile names -> outfile mappings
    As needed when file names do not map to a regular expression, specify
    file name, kit#, and surname in the corresponding csv file. There must be
    a kit# and name, but they can be generic strings like '000000'. If kit#
    and name do not form a unique tuple, there's a namespace clash, and a
    warning is issued.


Filename mapping spreadsheet (default name filemap.csv):

One mapping per line in csv. Order is unimportant. md5sum is optional. It is
not necessary to list files in the spreadsheet if they parse normally. It is
not necessary to list files with identical data more than once, even if the
filename differs, as long as the md5sum is correct.

Notes:
1. Known bad files can be mapped to None, None
2. There's a warning if multiple files map to the same name, but they have
   different data. Best to address all warnings before a production run
3. It's innocuous to have multiple file names, each with the same data,
   mapping to the same file name
4. It's possible and undesirable to create duplicate data if you map
   identical data to two different file names; this is reported in --dryrun
5. Duplicate lines may be left in the spreadsheet for clarity or may be pared
   down to a single line with the correct md5sum. Either way is functional.
   All identical files should be mapped to the same kit,name tuple.
6. Disable a mapping without removing it: put # in front of file name; then
   that line in the spreadsheet will not be applied
7. Every effort is made to have surnames and kit numbers exactly match other
   well-known places the data is published (esp., "The Big Tree", ytree.net)
8. When editing and importing/exporting in spreadsheet software, a string such
   as '000000' might get turned into a number (0). There are usually
   import/export options to ensure the column is text and avoid this

Column: filename
  The name of the original .zip. Applies to files that exactly match. Can be a
  blank value; then only matching via the md5sum applies (see
  below). filename+md5sum need form a unique value.

Column: kitid
  Typically the kit number. Can also be any other identifier that taken with
  the surname forms a unique value within the set of files.

Column: surname
  Surname or other descriptor such as Sweden or Unknown

Column: md5sum
  This does not need to be provided; however, it can help prevent ending up
  with some duplicates and inconsistently-named files when extracting data that
  contains more than one input zip (different names but same data). If a
  matching md5sum is found here, the corresponding naming in the dictionary is
  used. Beyond handling a file that has been renamed several times, another
  specific condition that is addressed can be explained by an example. Suppose
  the original .zip file was named "a.zip". A well-intentioned file maintainer
  came along and renamed this file to a more descriptive "fred-1234.zip."
  Further, Fredd, the file owner uploaded a second version of the same file
  named "Fredd-1234.zip." Now, someone trying to extract these files has both
  fred-1234.zip and Fredd-1234.zip but does not have a.zip. If a.zip is the
  only one of these in the dictionary above, upon extraction, the resulting
  file name is indeterminate. However, if the md5sum appears in the list,
  the algorithm will choose the dictionary entry, if it exists.  In this way,
  it doesn't matter which of the matching files you have, and only one of the
  three need be in the dictionary.

Other columns are ignored.
"""


import csv, os, sys, shutil, random
import re
import argparse
from time import time
from collections import defaultdict

# default - how verbose - the higher the number, the more chatty
DEBUG = 0

# process command line
parser = argparse.ArgumentParser()
parser.add_argument('-z', '--zipdir', required=True, action='append')
parser.add_argument('-u', '--unzipdir', required=True)
parser.add_argument('-d', '--dryrun', action='store_true')
parser.add_argument('-r', '--rename', action='store_true')
parser.add_argument('-v', '--verbose', action='count')
parser.add_argument('-k', '--keep', action='store_true')
parser.add_argument('-m', '--mapcsv', default='filemap.csv')
parser.add_argument('-s', '--subset', default='1:1')
namespace = parser.parse_args(sys.argv[1:])
command_args = {k:v for k,v in vars(namespace).items() if v}
zip_dirs = command_args['zipdir']
unzip_dir = command_args['unzipdir']
keep_files = vars(namespace)['keep']
verbose = vars(namespace)['verbose']
mapcsv = vars(namespace)['mapcsv']
keepbin = [s for s in vars(namespace)['subset'].split(':')]
keepbin = [int(keepbin[0]), [int(s) for s in keepbin[1].split(',')]]
if keepbin[1][0] > keepbin[0]:
    raise ValueError('check subset arg')

if not verbose:
    verbose = DEBUG

if vars(namespace)['dryrun']:
    if verbose < 1:
        verbose = 1

# trace (print) output if it exceeds the noise level
def trace (level, msg):
    if level <= verbose:
        print(msg)
        sys.stdout.flush()

# Always-apply rules for names. Surnames can be adjusted here. The only things
# that should be done in this routine are universal rules. For example, if,
# instead of "&#39;" in the surname, you always want apostrophe, you can do it
# here and keep some things out of the filemap.
def name_preference(name):
    newname = name.replace('&#39;', "'")
    # retain capitalization unless it's all upper or all lower case
    if newname.upper() == newname or newname.lower() == newname:
        newname = newname.title()
    return newname

# md5 hash of a given file name
def md5(fn):
    import hashlib
    md5hash = hashlib.md5()
    with open(fn, 'rb') as f:
        md5hash.update(f.read())
    return md5hash.hexdigest()

# hash function for subsetting - return True if this file is in the "kept" set
def keepfile(md5sum):
    m = int(md5sum,16)
    nbins = keepbin[0]
    if nbins <= 1:
        return True
    # special case - return a random set
    if keepbin[1][0] < 1:
        return (random.randint(1,keepbin[0]) == 1)
    # simple mod function is good enough hash here
    if (m % nbins + 1) in keepbin[1]:
        return True
    return False

# files that defy parsing the name are hand-crafted in this mapping
rename_dict = {}
with open (mapcsv) as csvfile:
    maps = csv.DictReader(csvfile)
    for row in maps:
        try:
            k = row['kitid'].strip()
            n = row['surname'].strip()
            m = row['md5sum'].strip()
            f = row['filename'].strip()
        except:
            trace(1, 'Warning: something is wrong with {}'.format(mapcsv))
        if f.startswith('#'):
            continue
        if f:
            rename_dict[f] = (k,n)
        if m:
            rename_dict[m] = (k,n)


trace(1, 'Verbosity: %d' % verbose)
if keep_files:
    trace(1, 'Keeping previously-extracted files - be sure they are OK!')

# create new directory for output
def setup_dirs (unzip_dir):
    if vars(namespace)['dryrun']:
        trace(1, 'Dry run: not making or purging {}'.format(unzip_dir))
        return
    if not keep_files:
        shutil.rmtree(unzip_dir, ignore_errors=True)
    try:
        os.makedirs(unzip_dir)
    except FileExistsError:
        pass
    except:
        print('ERROR: could not make directory', unzip_dir)

# extract zip files
# messy problem - messy solution - kit names are not consistent
def extract_zips(unzip_dir, files):

    # Handles de-duping; if there are multiple files with the same data,
    # we will use only one of the equivalent set. If there is an entry in
    # the filename dictionary, that takes precedence.
    md5sums = [md5(fn) for fn in files]
    mdict = defaultdict(list)
    for md,fn in zip(md5sums, files):
        mdict[md].append(fn)

    # Try to parse out at least the kit number by trying a series of regular
    # expressions. Adding regular expressions at the end of this list is safer
    # than at the beginning. Order is important - rules at top are matched
    # first.

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
    sepp = r'[\-\s\._]+'
    sept = r'[\-\s\._]{3}'
    bigy = r'(?:big' +seps+ r'y(?:data)?|ydna)'
    rslt = r'(?:results|data|rawdata|vcfdata|raw data|csvexport|raw_data|raw|bigyrawdata)'
    name = r'((?:'+nam2+seps+'){1,3})'
    kit = r'(?:(?:kit|ftdna)?[ #]?)?([enhb1-9][0-9]{3,6})'
    rzip = r'zip(?:.zip)?'
    snps = r'(?:[\-_]{,3}(?:r[\-_]?)?(?:cts\d{3,6}|fcg\d{4,5}|fgc\d{3,5}x?|p312|z\d{3,5}|df\d{2,3}x?|l\d{2,3}x?|u152|rs\d{4}|d27|sry\d{4}|m222|l\d{4}|s\d{4,6}|mc14|a\d{3,5}|zz\d{2}|zp\d{2}|z\d{2,3}|s\d{3}|pf\d{4}|by\d{3,5}|u106|l2|y\d{4,5}|yp\d{4,5})){1,3}'
    plac = r'(?:Germany|England|UnknownOrigin|Sweden|France|United_?Kingdom|Scotland|Ireland|Netherlands|Europe|Luxembour?g|Wales|Poland|Italy|CzechRepublic|Russia|Puerto-Rico|Switzerland|Algeria|Denmark|Slovakia|US|USA)'
    name_re = [
        #0 e.g. bigy-Treece-N4826.zip
        (re.compile(ws+sepp.join([bigy,name,kit,rzip]), re.I), 'name', 'kit'),
        #1 e.g. bigy-N4826-Treece.zip
        (re.compile(ws +sepp.join([bigy,kit,name,rzip]), re.I), 'kit', 'name'),
        #2 e.g. N4826-bigy-Treece.zip
        (re.compile(ws +sepp.join([kit,bigy,name,rzip]), re.I), 'kit', 'name'),
        #3 e.g. Treece - N4826 - bigy.zip
        (re.compile(ws+name  +sept+kit +sept +bigy +sep +r'?\.zip', re.I), 'name', 'kit'),
        #4 e.g. Treece N4826 bigy.zip
        (re.compile(ws +sepp.join([name,kit,bigy,rzip]), re.I), 'name', 'kit'),
        #5 e.g. Treece N4826 bigy results 20140808.zip
        (re.compile(ws +sepp.join([name,kit,bigy,rslt,ndate,rzip]), re.I), 'name', 'kit'),
        #6 e.g. bigy-Treece-N4826-FGC1233.zip
        (re.compile(ws +sepp.join([bigy,name,kit,snps,rzip]), re.I), 'name', 'kit'),
        #7 e.g. FGC1234-N4826-Treece-England-bigy-rawdata-20140708.zip
        (re.compile(ws +sepp.join([snps,kit,name,plac,bigy,rslt,ndate,rzip]), re.I), 'kit', 'name'),
        #8 e.g. FGC1234-N4826-Treece-bigy-rawdata-20140708.zip
        (re.compile(ws +sepp.join([snps,kit,name,bigy,rslt,ndate,rzip]), re.I), 'kit', 'name'),
        #9 e.g. FGC1234-Treece-N4826-bigy-rawdata-20140708.zip
        (re.compile(ws +sepp.join([snps,name,kit,bigy,rslt,ndate,rzip]), re.I), 'name', 'kit'),
        #10 e.g. FGC1234-N4826-Treece-bigy-rawdata.zip
        (re.compile(ws +sepp.join([snps,kit,name,bigy,rslt,rzip]), re.I), 'kit', 'name'),
        #11 e.g. FGC1234-N4826-Treece-bigy-rawdata.zip
        (re.compile(ws +sepp.join([snps,kit,name,bigy,rzip]), re.I), 'kit', 'name'),
        #12 e.g. FGC1234-N4826-Treece-England-bigy-rawdata.zip
        (re.compile(ws +sepp.join([snps,kit,name,plac,rslt,ndate,rzip]), re.I), 'kit', 'name'),
        #13 e.g. N4826_Treece_US_BigY_RawData_2018-01-03.zip
        (re.compile(ws +sepp.join([kit,name,plac,bigy,rslt,ndate,rzip]), re.I), 'kit', 'name'),
        #14 e.g. FGC1234-N4826-Treece-bigy-rawdata.zip
        (re.compile(ws +sepp.join([snps,kit,name,rslt,ndate,rzip]), re.I), 'kit', 'name'),
        #15 e.g. FGC1234-N4826-Treece-bigy-20140708.zip
        (re.compile(ws +sepp.join([snps,kit,name,bigy,ndate,rzip]), re.I), 'kit', 'name'),
        #16 e.g. FGC1234-N4826-Treece-20140708.zip
        (re.compile(ws +sepp.join([snps,kit,name,ndate,rzip]), re.I), 'kit', 'name'),
        #17 e.g. N4826-Treece-bigy-data-20140708.zip
        (re.compile(ws +sepp.join([kit,name,bigy,rslt,ndate,rzip]), re.I), 'kit', 'name'),
        #18 e.g. N4826-bigy-Treece-20140708.zip
        (re.compile(ws +sepp.join([kit,bigy,name,ndate,rzip]), re.I), 'kit', 'name'),
        #19 e.g. FGC1234-Treece-N4826.zip
        (re.compile(ws +sepp.join([snps,name,kit,rzip]), re.I), 'name', 'kit'),
        #20 e.g. FGC1234-Treece-N4826-bigy-rawdata.zip
        (re.compile(ws +sepp.join([snps,name,kit,bigy,rslt,rzip]), re.I), 'name', 'kit'),
        #21 e.g. bigy-Lindström-548872.zip
        (re.compile(ws+sepp.join([bigy,cname,kit,rzip]), re.I), 'name', 'kit'),
        ]

    trace (2, 'File names mapped, according to which regular expression:')
    # track counts - only for diagnostics
    cnt = defaultdict(int)
    # list of non-matching files
    nomatch=[]
    # all of the file names we could parse
    fname_dict = {}

    bed_re = re.compile(r'(\b(?:\w*[^_/])?regions(?:\[\d\])?\.bed)')
    vcf_re = re.compile(r'(\b(?:\w*[^_/])?variants(?:\[\d\])?\.vcf)')
    zip_re = re.compile(r'(\b(?:\w*[^_/])?bigy.*\.zip)')
    mgroup = []
    for md,fnames in mdict.items():
        if md in rename_dict:
            kkit, nname = rename_dict[md]
            pathname = fnames[0]
            fname = os.path.split(pathname)[-1]
            rule = 'dm'
        else:
            for pathname in fnames:
                fname = os.path.split(pathname)[-1]
                if fname in rename_dict:
                    kkit, nname = rename_dict[fname]
                    rule = 'df'
                    break
            else:
                kkit = nname = None
        if not keepfile(md):
            trace(2, '{} skipped because of subsetting'.format(fname))
            continue

        if kkit:
            if kkit in ('None', '') and nname in ('None', ''):
                kkit = nname = None
                trace(2, '{} skipped because of entry in mappings dictionary'.format(fname))
            else:
                trace(2, '{3:>2} {0:<50s} {1:<15s} {2:<10s}'.format(fname, nname, kkit, rule))
                cnt[rule] += 1
            fname_dict[pathname] = kkit, nname
        else:
            pathname = fnames[0] # only use one of the equivalent set
            if os.path.splitext(pathname)[-1].lower() != '.zip':
                trace (2, 'Found foreigner hanging out in zip directory: {0}'.format(fname))
                continue
            d = {}
            fname = os.path.split(pathname)[-1]
            for ii, (r,k1,k2) in enumerate(name_re):
                s = r.search(fname)
                if s:
                    d[k1] = s.groups()[0]
                    if k2:
                        d[k2] = s.groups()[1]
                    else:
                        d['name'] = 'Unknown'
                    d['name'] = name_preference(d['name'])
                    try:
                        trace (2, '{3:>2} {0:<50s} {1:<15s} {2:<10s}'.format(fname,
                                                   d['name'], d['kit'], ii))
                        cnt[ii] += 1
                        fname_dict[pathname] = d['kit'], d['name']
                    except:
                        trace (0, 'FAILURE on filename:', fname)
                    break
            else:
                nomatch.append(pathname)
        if len(fnames) > 1:
            for eq in [p for p in fnames if p!=pathname]:
                trace(2, '  -same: {}'.format(os.path.split(eq)[-1]))

    trace (2, 'Number of filenames not matched: {0}'.format(len(nomatch)))
    trace (2, 'Which expressions were matched:')
    def keyfunc(v):
        return '{0!s:0>2}'.format(v)
    for nn in sorted(cnt,key=keyfunc):
        trace (2, '{0:>2}: {1:>4}'.format(nn,cnt[nn]))

    if nomatch:
        trace (1, 'Files that did not match:')
        for ll in nomatch:
            trace (1, ll.strip())
    else:
        trace (1, 'All files matched a rule')

    zipcount = 0

    # keep track of what needs to be cleaned up
    emptydirs = []

    import zipfile

    for fname in fname_dict:
        kitnumber, kitname = fname_dict[fname]
        if kitnumber == None:
            continue
        if keep_files:
            vcffile = os.path.join(unzip_dir, '%s-%s.vcf' % (kitname, kitnumber))
            bedfile = os.path.join(unzip_dir, '%s-%s.bed' % (kitname, kitnumber))
            # no checking to see if the contents are good, but that's what was asked for
            if os.path.isfile(vcffile) and os.path.isfile(bedfile):
                trace(1, '%s-%s already exists and keep flag - skipping' % (kitname, kitnumber))
                continue
        try:
            zf = zipfile.ZipFile(fname)
        except:
            trace (0, 'WARN: not a zip file: {} from {}'.format(fname, os.getcwd()))
            continue
        listfiles = zf.namelist()
        bedfile = vcffile = zipfname = None
        for ff in listfiles:
            dirname, basename = os.path.split(ff)
            if bed_re.search(basename):
                bedfile = ff
            elif vcf_re.search(basename):
                vcffile = ff
            elif zip_re.search(basename):
                zipfname = ff
            if dirname and (dirname not in emptydirs):
                emptydirs.append(dirname)
        if (not bedfile) or (not vcffile):
            if not zipfname:
                trace(0, 'WARN: missing data in '+fname)
                continue
            else:
                try:
                    zf.extractall(unzip_dir, [zipfname,])
                    emptydirs.append(os.path.join(unzip_dir, zipfname))
                    bedfile = 'regions.bed'
                    vcffile = 'variants.vcf'
                    zf = zipfile.ZipFile(os.path.join(unzip_dir, zipfname))
                except:
                    trace(0, 'WARN: missing data in '+fname)
                    continue
        try:
            zf.extractall(unzip_dir, [bedfile, vcffile])
        except RuntimeError:
            trace(0, 'WARN: {} would not extract - encrypted?'.format(base))
        base = '%s-%s' % (kitname, kitnumber)

        fpath = os.path.join(unzip_dir, '%s')
        trace (3, fpath % base)
        if vars(namespace)['rename']:
            try:
                os.link(fname, (fpath % base)+'.zip')
                trace(1, 'ln {} {}.zip'.format(fname,(fpath % base)))
            except:
                shutil.copy2(fname, (fpath % base)+'.zip')
                trace(1, 'cp -p {} {}.zip'.format(fname,(fpath % base)))
            emptydirs.append(fpath % bedfile)
            emptydirs.append(fpath % vcffile)
        else:
            try:
                os.rename(fpath % bedfile, (fpath % base)+'.bed')
                os.rename(fpath % vcffile, (fpath % base)+'.vcf')
            except:
                trace(0, 'WARN: could not identify VCF and/or BED file for '+base)
        zipcount += 1

    trace (0, '%d new files extracted' % zipcount)

    # clean up any empty dirs unzip created
    if emptydirs:
        trace (3, 'Trying to remove droppings:')
        for dir in emptydirs:
            if os.path.isfile(dir):
                os.unlink(dir)
        for dir in emptydirs:
            if os.path.isfile(dir):
                continue
            try:
                dp = os.path.join(unzip_dir, dir)
                os.removedirs(dp)
                trace (3, '  {0}'.format(dp))
            except FileNotFoundError:
                pass
            except:
                trace (3, '  W! could not remove {0}'.format(dp))
                pass

    # list of file names we unzipped
    files = os.listdir(unzip_dir)
    return files

# some sanity checks for the files just extracted
def check_extraction(dirname):
    import stat, os, re, string

    trace(1, 'Checking names and checking for duplicate content...')
    files = os.listdir(dirname)
    os.chdir(dirname)
    fns = {}
    for f in files:
        if not (f.endswith('.bed') or f.endswith('.vcf')):
            trace(1, 'WARN badly named file {} found'.format(f))
            continue
        if f.endswith('.bed'):
            md5sum = md5(f)
            if md5sum in fns:
                trace(1, 'WARN: duplicate files: {},{}'.format(f,fns[md5sum]))

    trace(1, 'Checking for missing BED or VCF...')
    for f in files:
        bn,ext = os.path.splitext(f)
        if not os.path.isfile(bn + '.vcf') or not os.path.isfile(bn + '.bed'):
            trace(1, 'WARN: {} is missing a vcf or bed file'.format(f))
    trace(1, 'Checking for minimum file size and permissions...')
    for f in files:
        finfo = os.stat(f)
        if finfo.st_size < 100000 and f.endswith('.bed'):
            trace(1, 'WARN: {} is below minimum size for a BED'.format(f))
        elif finfo.st_size < 1000000 and f.endswith('.vcf'):
            trace(1, 'WARN: {} appears to be too small for a VCF'.format(f))
        if not (finfo.st_mode & stat.S_IRUSR):
            trace(1, 'WARN: {} is not readable.'.format(f))
    trace(1, 'Checking file names...')
    for f in files:
        if not f.isprintable():
            trace(1, 'WARN: {} is not printable.'.format(f))
        if re.search('['+string.whitespace+']', f):
            trace(1, 'WARN: {} has whitespace.'.format(f))


def main():
    import tempfile
    global unzip_dir
    # collect run time statistics
    trace(1,'Running the unpack-zip script...')
    T0 = time()
    files = []
    for zip_dir in zip_dirs:
        fl = [os.path.join(zip_dir,fn) for fn in os.listdir(zip_dir)]
        for fn in fl:
            if os.path.isfile(fn):
                files.append(fn)
    with tempfile.TemporaryDirectory() as tmpdirname:
        setup_dirs(unzip_dir)
        if vars(namespace)['dryrun']:
            unzip_dir = tmpdirname

        fnames = extract_zips(unzip_dir, files)
        trace (2, 'Number of files: {0}'.format(len(fnames)))
        trace (3, 'Files unpacked:')
        for ff in fnames:
            trace (3, ff)
        if vars(namespace)['dryrun'] and not vars(namespace)['rename']:
            trace(1, 'Running some sanity checks on extracted files...')
            check_extraction(tmpdirname)
    T1 = time()
    trace(1, '...complete after %f seconds' % (T1-T0))
    sys.exit(0)



main()
