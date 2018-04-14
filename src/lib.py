#!/usr/bin/env python3
# coding: utf-8
#
# Copyright (c) 2018 The Authors
#
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
#
# Usage:
#   import as a library; or
#   run script as a command with no args to execute some tests
#
# Environment:
#   REDUX_PATH must be set to source directory
#   config.yaml is read for configuration settings
#

import os, yaml, shutil, re, csv, zipfile, subprocess
from db import DB
import time
import sys
import hashlib
import requests, json
import urllib, time


# read the config file
sys.path.insert(0, os.environ['REDUX_PATH'])
REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
config = yaml.load(open(REDUX_CONF))

# add bin directory to path for executing utilities
# sys.path.insert(0, config['REDUX_BIN'])


# Class: DisplayMsg
# Purpose: mix-in used by Trace to display or log messages
class DisplayMsg:
    'class that provides methods to send messages to a file or a stream'
    def display (self, prefix='', msg='', stream=sys.stdout):
        'send a message to an output stream or file'
        print(prefix+msg, file=stream)
    def flush (self, stream=sys.stdout):
        'flush the output stream'
        stream.flush()

# Class: Trace
# Purpose: provides methods for debugging and status messages
# Initionalization variables:
#   prefix, optional text string prepended to every message
#   stream, optional output stream or file object
#   level, optional verbosity threshhold
#   autoflush, optional flag causes automatic flushing on each message
# Examples:
#   trace = Trace(1, stream=open('log.out', 'w'), autoflush=True)
#   trace(0, 'important message')
#   savetrace = trace
#   trace = Trace(99, stream=open('log.out', 'w'), prefix='dbg myfunc: ')
#   trace(10, 'local debugging message')
#   trace = savetrace
class Trace(DisplayMsg):
    'class that provides a logging facility: set level, then call trace method'
    def __init__ (self, level=0, prefix='', stream=sys.stderr, autoflush=False):
        self.level = level
        self.prefix = prefix
        self.stream = stream
        self.autoflush = autoflush
    def trace (self, level, msg):
        'trace(level,message): output message if level <= verbosity'
        if level <= self.level:
            self.display (self.prefix, msg, self.stream)
            if self.autoflush:
                self.flush()
    def __call__ (self, level, msg):
        self.trace (level, msg)

trace = Trace(config['verbosity'])

# Procedure: md5
# Purpose: return a md5 hash of a given object as a string signature
def md5(obj):
    md5hash = hashlib.md5()
    md5hash.update(str(obj).encode('utf-8'))
    return md5hash.hexdigest()

# Procedure: data_path
# Purpose: return a path to a file or directory in the configured data dir
def data_path(fname):
    return os.path.join(config['REDUX_DATA'], fname)

# Procedure: extract_zipdir
# Purpose:
#   call external utility to unpack zip files
#   this is driven from the list of DNA kits of interest
# Info:
#   this procedure is incomplete
#   the future intent is extract .vcf and .bed from .zip for debugging
#   currently, the extracted files are not needed; they do not land on disk
def extract_zipdir():
    # work in progress Jef
    # fixme - standalone utility - not needed for analysis
    # fixme - interop with API from DW
    # fixme - don't purge unzip dir if requested (-k flag)
    # fixme - regular expressions for the scripts need adjustment
    # callout to unpack-zip-files.py, utility from v1, in bin directory
    trace(10,'   Running the unpack-zip script...')
    io = subprocess.run(['unpack-zip-files.py', '-z', data_path(config['zip_dir']),
                             '-u', data_path(config['unzip_dir']),
                             '-m', data_path('filemap.csv')],
                            stdout=subprocess.PIPE)
    output = io.stdout.decode('utf-8') # ignore this output?
    fnames = os.listdir(data_path(config['unzip_dir']))
    trace (10, '   Number of files: {0}'.format(len(fnames)))
    trace (40, '   Files unpacked:')
    for ff in fnames:
        trace (40, ff)

# Procedure: getVCFvariants
# Input: an opened VCF file object
# Purpose: parse the vcf file for data to store in the database
def getVCFvariants(FILE):
    from subprocess import Popen, PIPE, STDOUT
    cmd = os.path.join(config['REDUX_ENV'], 'getVCFvariants.sh')
    p = Popen([cmd, '-'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    p_stdout = p.communicate(input=FILE.read())[0]
    output = p_stdout.decode('utf-8')
    # trace(1, 'command output: {}'.format(output))
    return output

# Procedure: readHg19Vcf
# Input: a VCF file name
# Purpose: read a hg19 FTDNA .zip file
# Returns: a dict of position -> mutation mappings
# Info:
#   This routine is not currently used because we focus on hg38. Modified from
#   Harald's analyzeVCF, this version returns every mutation with its derived
#   value, regardless of whether it was ancestral or not.
def readHg19Vcf(file):
    with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
        trace (30, "   Extracting VCF: %s" % vcffile)
        result = {}
        for line in vcffile:
            fields = line.split()
            if (fields[0] == 'chrY' and int(fields[1]) > 0 and fields[3] != '.'
                    and fields[4] != '.'):
                result[fields[1]] = [int(fields[1]), str(fields[3]),
                                         str(fields[4])]
        return result

# Procedure: populate_age
# Input: a database object
# Purpose: populate agebed table from age.bed (a set of ranges in a file)
# Returns: nothing
# Info:
#   this procedure discards what may be already in the agebed table and
#   replaces it with the ranges defined in age.bed, a text file
def populate_age(dbo):
    trace(1, 'populate agebed table')
    with open('age.bed') as bedfile:
        cf = csv.reader(bedfile, delimiter=' ')
        ranges = []
        for row in cf:
            try:
                ranges.append((row[0], row[1], row[2]))
            except:
                trace(0,'failed on row of age.bed:{}'.format(row))
        dc = dbo.cursor()
        dc.execute('delete from agebed')
        dc.execute('drop table if exists tmpt')
        dc.execute('create temporary table tmpt(a,b,c)')
        dc.executemany('insert into tmpt values(?,?,?)', ranges)
        dc.execute('''insert or ignore into bedranges(minaddr,maxaddr)
                           select b,c from tmpt''')
        dc.execute('''insert into agebed(bID)
                           select id from bedranges
                           inner join tmpt t on t.b=minaddr and t.c=maxaddr''')
        dc.execute('drop table tmpt')
    return

# Procedure: populate_refpos
# Purpose: fill up the refpos (reference-positive) table
# Input: a database object
# Returns: nothing
# Info:
#   populate reference positives from the SNPs listed in refpos.txt, a text
#   file
def populate_refpos(dbo):
    trace(1, 'populate refpos table')
    with open('refpos.txt') as refposfile:
        cf = csv.reader(refposfile)
        snps = []
        for row in cf:
            if row[0].startswith('#'):
                continue
            try:
                snps.append(row[0])
            except:
                trace(0, 'failed on row of refpos.txt:{}'.format(row))
        dc = dbo.cursor()
        dc.execute('delete from refpos')
        for snpname in snps:
           dc.execute('''insert or ignore into refpos
                       select vid from snpnames where snpname=?''',
                       (snpname,))

        # make sure there's a swapped variant
        dc.execute('''insert or ignore into variants(pos,anc,der,buildid)
                      select pos,der,anc,buildid from variants v
                      inner join refpos r on r.vid=v.id''')
    return

# Procedure: populate_analysis_kits
# Purpose: fill up the list of kits to be used for analysis
# Input: a database object
# Returns: nothing
# Info:
#   populate analysis kits from those listed in kits.txt, a text file
def populate_analysis_kits(dbo):
    trace(1, 'populate analysis kit list table')
    with open('kits.txt') as kitsfile:
        cf = csv.reader(kitsfile)
        kitids = []
        for row in cf:
            if row[0].startswith('#'):
                continue
            try:
                kitids.append(row[0])
            except:
                trace(0, 'failed on row of kits.txt:{}'.format(row))
        dc = dbo.cursor()
        dc.execute('delete from analysis_kits')
        for kit in kitids:
           dc.execute('''insert into analysis_kits
                       select d.DNAID from dataset d
                       where d.kitID like ?''',
                       (kit,))
    return


# Procedure: populate_STRS
# Purpose: populate a table of STR definitions
# Input:
#   a database object
#   optional ordering (default order is same as FTDNA)
# Info:
#   This table is not yet used; we're setting up a framework for using STR data
#   to incorporate into the models.  Ordering is optional. If not given, table
#   stores SNP names in FTDNA Y111 order.  Otherwise, It's a 111-number vector.
def populate_STRs(dbo, ordering=None):
    strdefs = (
        'DYS393', 'DYS390', 'DYS19', 'DYS391', 'DYS385a', 'DYS385b',
        'DYS426', 'DYS388', 'DYS439', 'DYS389i', 'DYS392', 'DYS389ii',
        'DYS458', 'DYS459a', 'DYS459b', 'DYS455', 'DYS454', 'DYS447',
        'DYS437', 'DYS448', 'DYS449', 'DYS464a', 'DYS464b', 'DYS464c',
        'DYS464d', 'DYS460', 'YH4', 'YCAIIa', 'YCAIIb', 'DYS456', 'DYS607',
        'DYS576', 'DYS570', 'CDYa', 'CDYb', 'DYS442', 'DYS438', 'DYS531',
        'DYS578', 'DYF395S1a', 'DYF395S1b', 'DYS590', 'DYS537', 'DYS641',
        'DYS472', 'DYF406S1', 'DYS511', 'DYS425', 'DYS413a', 'DYS413b',
        'DYS557', 'DYS594', 'DYS436', 'DYS490', 'DYS534', 'DYS450',
        'DYS444', 'DYS481', 'DYS520', 'DYS446', 'DYS617', 'DYS568',
        'DYS487', 'DYS572', 'DYS640', 'DYS492', 'DYS565', 'DYS710',
        'DYS485', 'DYS632', 'DYS495', 'DYS540', 'DYS714', 'DYS716',
        'DYS717', 'DYS505', 'DYS556', 'DYS549', 'DYS589', 'DYS522',
        'DYS494', 'DYS533', 'DYS636', 'DYS575', 'DYS638', 'DYS462',
        'DYS452', 'DYS445', 'YA10', 'DYS463', 'DYS441', 'Y1B07', 'DYS525',
        'DYS712', 'DYS593', 'DYS650', 'DYS532', 'DYS715', 'DYS504',
        'DYS513', 'DYS561', 'DYS552', 'DYS726', 'DYS635', 'DYS587',
        'DYS643', 'DYS497', 'DYS510', 'DYS434', 'DYS461', 'DYS435')
    dc = dbo.cursor()
    if ordering and (len(ordering)==len(strdefs)):
        for tup in zip(ordering,strdefs):
            dc.execute('insert or ignore into strs(ordering,strname) values(?,?)', tup)
    else:
        for tup in enumerate(strdefs):
            dc.execute('insert or ignore into strs(ordering,strname) values(?,?)', tup)

# Procedure: get_kits
# Purpose: pull information about the kits from the web api of haplogroup-r
# Input:
#   fromweb: if True, try to refresh from the web; else, prefer cached
# Returns:
#   js, a json object with all of the metadata records
# Info:
#   haplogroup-r provides an api to get metadata about kits stored there.  This
#   procedure calls that API and returns the JSON record.  To avoid making too
#   many repeated calls to the API when testing and developing, the json record
#   is cached on disk. When the API parameter is empty, satisfy the request
#   from the cached copy.
def get_kits (fromweb=True):
    # if pulling from the web, where to get the result
    API = 'http://haplogroup-r.org/api/v1/uploads.php'
    qry = 'format=json'
    # if re-using cached data from a previous web pull, where to find the file
    fname = data_path(os.path.join('cache','dataset.json'))
    if not fromweb:
        try:
            trace(1, 'reading kit info from {}'.format(fname))
            js = json.loads(open(fname).read())
            return js
        except:
            trace(0, 'no cached {} - trying web'.format(fname))

    try:
        trace(1, 'reading kit info from the web')
        url = '?'.join([API, qry])
        res = requests.get(url)
        js = res.json()
        open(fname,'w').write(json.dumps(js))
        return js
    except:
        trace(0, 'Failed to pull kit metadata from {}'.format(API))
        raise # fixme - what to do on error?

# Procedure: update_metadata
# Purpose: update the kit information in the database
# Input:
#   db, a database object
#   js, a json object that comes from the Haplogroup-R DW API
# Returns: nothing - only updates the database tables dataset, person, etc
# Info:
#   this procedure doesn't load any data; it just updates the metadata for the
#   available kits contained in the json record
def update_metadata(db, js):
    dc = db.cursor()
    blds = {'b38': 'hg38', 'b19': 'hg19', 'b37': 'hg19'}
    rows = [(
        # fields that go into dataset table directly
        jr['kitId'].strip(), jr['uploaded'], jr['dataFile'], jr['long'],
        jr['lat'], jr['otherInfo'], jr['origFileName'], jr['birthYear'],
        jr['approxHg'],
        # need inserting into their own tables
        jr['country'], jr['normalOrig'], jr['lab'], blds[jr['build']],
        jr['surname'], jr['testType'], jr['isNGS'])
        for jr in js
        ]
    trace(3, 'first row out of {}:{}'.format(len(rows),rows[0]))
    trace(1, '{} unique kit ids'.format(len(set([v['kitId'] for v in js]))))
    trace(1, '{} null surname'.format(len([v['surname'] for v in js if not v['surname']])))

    # populate the dependency tables
    # (testtype,isNGS) goes into testtypes
    tups = [y for y in set([(r[-2],r[-1]) for r in rows])]
    dc.executemany('insert or ignore into testtype(testNm,isNGS) values(?,?)', tups)

    # (surname+kitId+build) goes into person
    # fixme: we need a DNA-to-person mapping. This is a big kludge
    tups = [y for y in set([(r[-3],r[0],r[-4]) for r in rows])]
    dc.executemany('insert or ignore into person(surname,firstname,middlename) values(?,?,?)',
                          tups)

    for tbl,val,idx in (('country','country',-7),
                        ('surname', 'surname',-3),
                        ('origin', 'origin',-6),
                        ('lab', 'labNm',-5),
                        ('build', 'buildNm',-4)):
        tups = [(y,) for y in set([v[idx] for v in rows])]
        trace(3, 'first tuple to insert into {}: {}'.format(tbl, tups[0]))
        dc.executemany('insert or ignore into {}({}) values(?)'.format(tbl,val),
                              tups)

    # create temporary table, where columns correspond to values above
    dc.execute('''create temporary table tmpt(
           a TEXT, b TEXT, c TEXT, d TEXT,
           e TEXT, f TEXT, g TEXT, h TEXT,
           i TEXT,
           j TEXT, k TEXT, l TEXT, m TEXT,
           n TEXT, o TEXT, p TEXT)''')
    dc.executemany('''INSERT INTO tmpt(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', rows)
    dc.execute('''
        INSERT or ignore INTO dataset(kitId, importDt, fileNm, lng,
            lat, otherInfo, origFileNm, birthYr,
            approxHg,
            countryID, normalOrigID, labID, buildID, testTypeID, DNAID, surnameID)
        SELECT a, b, c, d,
            e, f, g, h,
            i,
            cn.id, oc.id, ln.id, bn.id, tt.id, pn.id, sn.id
        FROM tmpt
        INNER JOIN country cn ON
            tmpt.j = cn.country or (cn.country is NULL and tmpt.j is NULL)
        INNER JOIN origin oc ON
            tmpt.k = oc.origin or (oc.origin is NULL and tmpt.k is NULL)
        INNER JOIN lab ln ON tmpt.l = ln.labNm
        INNER JOIN build bn ON tmpt.m = bn.buildNm
        INNER JOIN surname sn ON
            tmpt.n = sn.surname or (sn.surname is NULL and tmpt.n is NULL)
        INNER JOIN person pn ON
            (tmpt.n = pn.surname or (tmpt.n is NULL and pn.surname is NULL)) AND
            (tmpt.a = pn.firstname or (tmpt.a is NULL and pn.firstname is NULL)) AND
            (tmpt.m = pn.middlename or (tmpt.m is NULL and pn.middlename is NULL))
        INNER JOIN testtype tt ON tmpt.o = tt.testNm and tmpt.p = tt.isNGS''')
    trace(1, '{} rows inserted into dataset'.format(dc.execute('select count(*) from dataset').fetchone()[0]))
    dc.execute('drop table tmpt')
    return

# Procedure: get_build_byname
# Purpose: get build identifier by its name; creates new entry if needed
# Input:
#   db, a database object
#   buildname, optional name of build, e.g. 'hg38'
# Info:
#   known aliases are reduced to one entry
def get_build_byname(db, buildname='hg38'):
    if buildname.lower().strip() in ('hg19', 'grch37', 'b19', 'b37'):
        buildname = 'hg19'
    elif buildname.lower().strip() in ('hg38', 'grch38', 'b38'):
        buildname = 'hg38'
    dc = db.cursor()
    dc.execute('select id from build where buildNm=?', (buildname,))
    bid = None
    for bid, in dc:
        continue
    if not bid:
        dc.execute('insert into build(buildNm) values (?)', (buildname,))
        bid = dc.lastrowid
    return bid

# Procedure: populate_contigs
# Purpose: load data into the contig table
# Input: a database object
def populate_contigs(db):
    trace(1,'populate contigs table')
    bid = get_build_byname(db, 'hg19')
    dc = db.cursor()
    dc.execute('insert into Contig(buildID,description,length) values(?,?,?)',
                      (bid, 'chrY', 59373566))
    bid = get_build_byname(db, 'hg38')
    dc.execute('insert into Contig(buildID,description,length) values(?,?,?)',
                      (bid, 'chrY', 57227415))
    return

# Procedure: populate_fileinfo
# Purpose: populate dataset information from Haplogroup-R data warehouse API
# Input:
#   db, a database object
#   fromweb, if True then request a pull of fresh data via the web
# Info:
#   This is all that needs to be called to update metadata about available kits
#   from haplogroup-r. In the normal case, it reaches out using the web API and
#   stores the latest data. It can also store the most-recent cached data
#   without calling the web API
def populate_fileinfo(dbo, fromweb=True):
    js = get_kits(fromweb)
    update_metadata(dbo, js)

# Procedure: updatesnps
# Purpose: update snp definitions
# Input:
#   db, a database instance
#   snp_reference, a csv DictReader instance with the SNP definitions
#   buildname, (optional, default 'hg38') the reference build name
# Returns: nothing - only updates the table
def updatesnps(db, snp_reference, buildname='hg38'):
    trace(1, 'update snpnames for {}'.format(buildname))
    bid = get_build_byname(db, buildname)
    dc = db.cursor()
    dc.execute('drop table if exists tmpt')
    dc.execute('create temporary table tmpt(a integer, b integer, c text, d text, e text, unique(a,b,c,d,e))')
    dc.executemany('INSERT OR IGNORE INTO tmpt(a,b,c,d,e) VALUES (?,?,?,?,?)',
             ((bid, rec['start'],
                   rec['allele_anc'].strip(), rec['allele_der'].strip(), rec['Name'])
                     for rec in snp_reference))
    dc.execute('''insert or ignore into alleles(allele) 
                        select distinct c from tmpt''')
    dc.execute('''insert or ignore into alleles(allele) 
                        select distinct d from tmpt''')
    dc.execute('''insert or ignore into variants(buildID, pos, anc, der)
                        select a, b, an.id, dr.id from tmpt
                        inner join alleles an on an.allele = c
                        inner join alleles dr on dr.allele = d''')
    dc.execute('''insert or ignore into snpnames(snpname,vID)
                        select t.e, v.id from tmpt t, variants v, alleles a, alleles d
                        where t.c=a.allele and v.anc=a.id and
                              t.d=d.allele and v.der=d.id and
                              v.pos = t.b and v.buildID=t.a''')

    dc.execute('drop table tmpt')
    return

# Procedure: get_SNPdefs_fromweb
# Purpose: pull SNP definitions from the web at ybrowse.org
# Input:
#   db, a database instance
#   maxage, maximum number of days since last update
#   url (optional), where to get the data file
# Info:
#   This should be called after the build table is populated. Uses urllib to
#   fetch the SNP definition file from the web if the latest version is older
#   than maxage. This keeps the SNP definitions reasonably up to date without
#   having to download the large files repeatedly. Refresh files if they are
#   older than maxage; do nothing if maxage < 0
def get_SNPdefs_fromweb(db, maxage, url='http://ybrowse.org/gbrowse2/gff'):
    UpdatedFlag = False
    if maxage < 0:
        return
    # convert to seconds
    maxage = 24*3600*maxage
    dc = db.cursor()
    for (build,) in dc.execute('select buildNm from build'):
        deltat = maxage + 1
        fbase = 'snps_{}.csv'.format(build)
        fname = data_path(os.path.join('cache', fbase))
        fget = os.path.join(url, fbase)
        try:
            if os.path.exists(fname):
                deltat = time.clock() - os.path.getmtime(fname)
            if deltat > maxage:
                trace (1, 'refresh: {}'.format(fbase))
                urllib.request.urlretrieve(fget, fname)
                deltat = time.clock() - os.path.getmtime(fname)
                UpdatedFlag = True
        except:
            pass
        if not os.path.exists(fname) or deltat > maxage:
            trace(0, 'failed to update {} from the web'.format(fname))
    return UpdatedFlag

# Procedure: populate_snps
# Purpose: populate SNP definitions in the database
# Input:
#   dbo, a database object
#   maxage (optional), maximum age of data files pulled from web
# Info:
#   refresh from web if we have is older than maxage (in days)
def populate_SNPs(dbo, maxage=config['max_snpdef_age']):
    # don't update if nothing changed
    if not (config['drop_tables'] or get_SNPdefs_fromweb(dbo, maxage=maxage)):
        return
    # update known snps for hg19 and hg38
    cachedir = os.path.join(config['REDUX_DATA'], 'cache')
    with open(os.path.join(cachedir, config['b37_snp_file'])) as snpfile:
        snp_reference = csv.DictReader(snpfile)
        updatesnps(dbo, snp_reference, 'hg19')
    with open(os.path.join(cachedir, config['b38_snp_file'])) as snpfile:
        snp_reference = csv.DictReader(snpfile)
        updatesnps(dbo, snp_reference, 'hg38')
    return


# Procedure: populate_from_BED_file
# Purpose: populate regions from a FTDNA BED file
# Input:
#   dbo, a database object
#   pid, a database person ID
#   fileobj, a file object from the open .bed file
def populate_from_BED_file(dbo, pid, fileobj):
    dc = dbo.cursor()
    ranges = []
    try:
        for line in fileobj:
            ychr, minr, maxr = line.split()
            ranges.append((pid, int(minr), int(maxr)))
    except:
        trace(0, 'FAILED on file at {}'.format(fileobj.readline()))
        return
    trace(500, '{} ranges for pID {}'.format(len(ranges), pid))

    dc.execute('drop table if exists tmpt')
    dc.execute('create temporary table tmpt(a,b,c)')
    dc.executemany('insert into tmpt values(?,?,?)', ranges)
    dc.execute('''insert or ignore into bedranges(minaddr,maxaddr)
                  select b,c from tmpt''')
    dc.execute('''insert into bed(pID, bID)
                  select t.a, br.id from bedranges br
                  inner join tmpt t on
                  t.b=br.minaddr and t.c=br.maxaddr''')
    dc.close()
    return

# Procedure: pack_call
# Purpose: pack call information into an integer
# Info:
#   this procedure exists to make vcfcalls table more compact
#   stores: pos, anc, der, passfail, BQ, MQ, nreads, passrate, gt
#   needs corresponding unpack_call
def pack_call(call_tup):
    maxshift = 33
    ngbits = 2
    nqbits = 7
    ncbits = 9
    npbits = 8
    # pack in a pass/fail flag
    bitfield = {'PASS': 1<<maxshift, 'FAIL': 0}[call_tup[3]]
    maxshift -= 1
    # pack in a genotype, [0,3]
    gtdict = {'0/0':0, '1/1':1, '0/2':2, '0/1':2, '1/2':3, '1/3':3, '2/2':3}
    bitfield |= (gtdict[call_tup[8]] << (maxshift - ngbits + 1))
    maxshift -= ngbits
    # pack in a unitless nqbits-number representing q1
    qnum = int(float(call_tup[4]) * 3)
    if qnum > (1<<nqbits) - 1:
        qnum = (1<<nqbits)-1
    bitfield |= (qnum << (maxshift - nqbits + 1))
    maxshift -= nqbits
    # pack in a unitless nqbits-number representing q2
    qnum = int(float(call_tup[5]) * 2)
    if qnum > (1<<nqbits) - 1:
        qnum = (1<<nqbits) - 1
    bitfield |= (qnum << (maxshift - nqbits + 1))
    maxshift -= nqbits
    # pack in a ncbits-bit number representing num reads
    nreads = int(call_tup[6])
    if nreads > (1<<ncbits) - 1:
        nreads = (1<<ncbits) - 1
    bitfield |= (nreads << (maxshift - ncbits + 1))
    maxshift -= nqbits
    # pack in a ncbits-bit number representing passrate
    passrate = int(float(call_tup[7]) * ((1<<npbits) - 1))
    if passrate > (1<<npbits) - 1:
        passrate = (1<<npbits) - 1
    bitfield |= passrate
    return bitfield

# Procedure: unpack_call
# Purpose: unpack callinfo (corresponds to pack_call)
# Info:
#   if pack_call changes, this needs to change
def unpack_call(bitfield):
    maxshift = 33
    ngbits = 2
    nqbits = 7
    ncbits = 9
    npbits = 8
    passrate = float(bitfield & (1<<npbits) - 1) / ((1<<npbits) - 1)
    calls = (bitfield>>npbits) & ((1<<ncbits) - 1)
    q2 = float((bitfield>>(npbits+ncbits)) & ((1<<nqbits) - 1)) / 2.
    q1 = float((bitfield>>(npbits+ncbits+nqbits)) & ((1<<nqbits) - 1)) /3.
    gt = (bitfield>>(npbits+ncbits+nqbits*2)) & ((1<<ngbits) - 1)
    if (1<<maxshift) & bitfield:
        passfail = True
    else:
        passfail = False
    return passfail, gt, q1, q2, calls, passrate

# Procedure: populate_from_VCF_file
# Purpose: populate calls, quality, and variants from a VCF file
# Input:
#   dbo, a database object
#   bid, a build ID
#   pid, a person ID
#   fileobj, a file object for the open VCF file
# Returns: nothing, only updates database tables
# Info:
#   depends on refpos table, which should already be populated
def populate_from_VCF_file(dbo, bid, pid, fileobj):
    dc = dbo.cursor()
    b = dc.execute('select buildNm from build where id=?', (bid,)).fetchone()[0]
    if b != 'hg38':
        trace(0, 'ERROR: currently unable to parse build {}'.format(b))
        return
    # parse output of getVCFvariants(fname)
    parsed = getVCFvariants(fileobj)
    #if pid == 1040:
    #    print(parsed)
    #    sys.exit()
    tups = []
    for line in parsed.splitlines():
        # pos, anc, der, passfail, q1, q2, nreads, passrate, gt
        tups.append(line.split())
    trace(5, 'parsed vcf: {}...'.format(tups[:3]))

    # filter down to the calls we want to store
    try:
        #passes = [t for t in tups if t[2] != '.' and (t[3] == 'PASS' or
        #          (int(t[6]) < 4 and float(t[7]) > .75))]
        passes = [t for t in tups]
    except ValueError:
        trace(0, 'parsing VCF failed on {}'.format(t))

    # save the distinct alleles - check performance
    alleles = set([x[1] for x in passes] + [x[2] for x in passes])
    dc.executemany('insert or ignore into alleles(allele) values(?)',
                       [(x,) for x in alleles])

    # save the call quality info - experimental
    call_info = [t + [pack_call(t)] for t in passes]

    # execute sql on results to save in vcfcalls
    dc.execute('drop table if exists tmpt')
    # bid, pos, ref, alt, pid, passfail, gt, packcall
    dc.execute('''create temporary table tmpt(a integer, b integer, c text,
                  d text, e integer, f text, g text, h integer)''')
    dc.executemany('insert into tmpt values(?,?,?,?,?,?,?,?)',
                    [[bid]+v[0:3]+[pid]+[v[3]]+v[-2:] for v in call_info])

    # don't insert clear reference calls unless they are refpos
    # FIXME - refpos test is not needed here because refpos variants do not
    # show up in the .vcf with derived = "." (but this line is not harmful)
    dc.execute('''delete from tmpt where g='0/0' and d='.' and f='PASS'
                  and b not in (select v.pos from variants v
                    inner join refpos r on r.vid=v.id)''')

    dc.execute('''insert or ignore into variants(buildID, pos, anc, der)
                  select a, b, an.id, dr.id from tmpt
                  inner join alleles an on an.allele = c
                  inner join alleles dr on dr.allele = d''')


    # fixme - performance?
    trace(4,'VCF update variants at {}'.format(time.clock()))

    dc.execute('''insert into vcfcalls (pid,vid,callinfo)
                  select e, v.id, h from tmpt
                  inner join alleles an on an.allele = c
                  inner join alleles dr on dr.allele = d
                  inner join variants v on v.buildID = a and v.pos = b
                  and v.anc=an.id and v.der=dr.id''')

    trace(4,'VCF load for {} done at {}'.format(pid, time.clock()))
    dc.execute('drop table tmpt')
    dc.close()

    trace(500, '{} vcf calls: {}...'.format(len(passes), passes[:5]))
    return

# Procedure: populate_from_zip_file
# Purpose: unpack any zip from FTDNA that has the bed file and vcf file
# Info:
#   INCOMPLETE
#   future: to permit import of .zip that isn't in Hap-R repository
def populate_from_zip_file(dbo, fname):
    # stub work in progress Jef
    # unpack bed and vcf to a temporary unzip dir
    # use the surname-kitID of unzip file as the kit name in the dataset table
    # create dataset entry if needed
    # pid = get_person_id
    # auto-determine build name for FTDNA to get build id (parsed from first lines)
    # bid = get_build_byname()
    # fname = unpacked vcf from unzip dir
    # populate_from_VCF_file
    # populate_from_BED_file
    return

# Procedure: populate_from_dataset
# Purpose:
#   unpack all zip files that can be handled from the HaplogroupR catalog.
# Input:
#   dbo, a database object
#   uses any zip files already downloaded into the data directory
# Info:
#   Uses previously-populated dataset table to find these files.
#   Assume we've already pulled the list from the H-R web API and we
#   already downloaded some zip files.
#
#   Loop over kit metadata in dataset table; if zip file exists locally:
#     if already loaded (calls exist in the db), skip this file
#     read the BED and VCF file from it in place, without landing on disk
#     store the VCF data in vcfcalls
#     store the BED ranges
#
# Future:
#   also download files?
def populate_from_dataset(dbo):
    trace(1, 'populate from dataset with kit limit {}'.format(config['kitlimit']))
    dc = dbo.cursor()
    bed_re = re.compile(r'(\b(?:\w*[^_/])?regions(?:\[\d\])?\.bed)')
    vcf_re = re.compile(r'(\b(?:\w*[^_/])?variants(?:\[\d\])?\.vcf)')
    # Query to get all known kits, with analysis_kits prioritized.
    # Prioritizing analysis_kits means they're loaded first, and we don't need
    # to load thousands of kits to get the ones we're interested in.
    dc.execute('''select fileNm,buildID,DNAID,1 from dataset
                  inner join analysis_kits on pID=DNAID
                        union all
                  select fileNm,buildID,DNAID,2 from dataset
                  order by 4''')
    allsets = list([(t[0],t[1],t[2]) for t in dc])
    pc = dbo.cursor()
    pl = pc.execute('select distinct pid from vcfcalls')
    pexists = [p[0] for p in pl]
    trace(5,'allsets: {}'.format(allsets[:config['kitlimit']]))
    nkits = 0

    # FIXME - hack - better index handling might be desirable.
    # Here we drop indexes while updating because inserts are faster.
    # Later, we re-create these indexes after all of the updates are done.
    trace(3, 'drop indexes at {}'.format(time.clock()))
    dc.execute('drop index if exists bedidx')
    dc.execute('drop index if exists vcfidx')
    dc.execute('drop index if exists vcfpidx')
    trace(3, 'done at {}'.format(time.clock()))

    for (fn,buildid,pid) in allsets:

        # begin work so we can roll back this kit if it fails
        trace(3, 'committing work')
        dbo.commit()

        # Loop over the kits we know about.
        # If there are already calls for this person, skip the load.
        if (not config['drop_tables']) and pid in pexists:
            trace(2, 'calls exist - skip {}'.format(fn[:50]))
            continue
        zipf = os.path.join(data_path('HaplogroupR'), fn)
        if not os.path.exists(zipf):
            trace(10, 'not present: {}'.format(zipf))
            continue
        try:
            # open the zip file and pull out the BED and VCF
            with zipfile.ZipFile(zipf) as zf:
                trace(1, '{}-{}'.format(nkits,zf.filename[:70]))
                listfiles = zf.namelist()
                bedfile = vcffile = None
                for ff in listfiles:
                    dirname, basename = os.path.split(ff)
                    if bed_re.search(basename):
                        bedfile = ff
                    elif vcf_re.search(basename):
                        vcffile = ff
                if (not bedfile) or (not vcffile):
                    trace(0, 'FAIL: missing data:{} (not loaded)'.format(zipf))
                    continue
                with zf.open(bedfile,'r') as bedf:
                    trace(3, 'populate from bed {}'.format(bedfile))
                    populate_from_BED_file(dbo, pid, bedf)
                with zf.open(vcffile,'r') as vcff:
                    trace(3, 'populate from vcf {}'.format(vcffile))
                    populate_from_VCF_file(dbo, buildid, pid, vcff)
            nkits += 1
        except:
            # something failed while loading this file - roll back changes
            trace(0, 'FAIL on file {} (not fully loaded)'.format(zipf))
            trace(0, 'roll-back this file load and continue')
            dbo.db.rollback()

        if nkits >= config['kitlimit']:
            break

    # re-create indexes we dropped above
    trace(3, 're-create indexes at {}'.format(time.clock()))
    dbo.commit()
    dc.execute('create index bedidx on bed(pID,bID)')
    dc.execute('create index vcfidx on vcfcalls(vID)')
    dc.execute('create index vcfpidx on vcfcalls(pID)')
    dc.close()
    trace(3, 'done at {}'.format(time.clock()))
    return

# Procedure: db_creation
# Purpose: a high-level routine that can be called for initial database
#   creation and table loads
def db_creation():
    db = DB(drop=config['drop_tables'])
    if config['drop_tables']:
        db.create_schema()
        populate_fileinfo(db, fromweb=config['use_web_api'])
        populate_STRs(db)
        populate_SNPs(db)
        populate_contigs(db)
        populate_age(db)
        populate_refpos(db)
        populate_analysis_kits(db)
    return db


# test framework
if __name__=='__main__':
    # unit test the Trace class
    t = Trace(1)
    t(1, 'test message should display to stdout')
    t(2, 'this message should not be seen')
    savet = t
    t = Trace(99, stream=open('log.out', 'w'), prefix='dbg myfunc: ')
    t.autoflush = True
    t(0, 'local debugging message should be in log.out')
    t(98, 'local debugging message should also be in log.out')
    t(100, 'this message should not be seen')
    t = savet
    t(0, 'another message that should be seen')
