#!/usr/bin/env python3
# coding: utf-8

# Copyright (c) 2018 The Authors

# Contributors: Jef Treece, Harald Alvestrand, Zak Jones, Iain McDonald
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# run this script as a command with no args to execute some tests
# it probably only works when run from the src directory


import os,yaml,shutil,glob,re,csv,zipfile,subprocess
from db import DB
from collections import defaultdict
from array_api import *
import time

# read the config file
# todo: this might need to some bootstrapping (run outside of src dir?)
config = yaml.load(open('config.yaml'))

# add src and bin directories to path
import sys
sys.path.insert(0, config['REDUX_PATH'])
sys.path.insert(0, config['REDUX_BIN'])


# routines - debug/diagnostic output
# fixme - there are too many levels of verbosity - probably should be a
# bitmap/flags
# output to stderr by default - use 2> redirection in bash to capture
# this sends all zero-level messages to stdout
def trace (level, msg, stream=sys.stderr):
    if level <= config['verbosity']:
        if level == 0:
            print(msg)
        else:
            print(msg, file=stream)
            stream.flush()

# return a path to a file or directory in the configured data dir
def data_path(fname):
    return os.path.join(config['REDUX_DATA'], fname)

def refresh_dir(DIR,cleanFlag=False):
    DIR = os.path.join(config['REDUX_ENV'], DIR)
    #print DIR
    if (os.path.isdir(DIR)):
        files = glob.glob(DIR+'/*')
        if cleanFlag:
            for f in files:
                os.remove(f)
    else:
        os.makedirs(DIR)

# update the mtime of a file without changing file contents
def touch_file(FILE):
    if os.path.exists(FILE):
        open(FILE,'a').close()

def cmd_exists(CMD):
    return any(os.access(os.path.join(path, CMD), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))


# various scratch/data directory setup under "data"
def setup_dirs():
    shutil.rmtree(data_path(config['unzip_dir']),ignore_errors=True)
    os.makedirs(data_path(config['unzip_dir']))

# call external utility to unpack all files in the zip directory
# return list of files that were unzipped
# this procedure is probably obsolete - remove in the future
def extract_zipdir():
    # work in progress Jef
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

# unpack a single zip file
def unpack_zipfile():
    # work in progress Jef
    return

    
# cache previous run's results
# fixme this is incorrect, and we need to figure out what needs to be saved
def go_backup():

    trace(0,"** performing backup.")
    trace(0,"** (msg) CREATING BACKUP COPIES OF EXISTING FILES...")
    
    # autobackup dir
    refresh_dir('autobackup')
    for FILE_PATTERN in config['backup_files'].split():
        for FILE in glob.glob(FILE_PATTERN):
            shutil.copy(FILE,'autobackup')
    
    if config['make_report']:
        #print "MAKING REPORT..."
        delete_file('report.csv')
    
    trace(0,"** + backup done.")
    


# utility procuedure to drop/create the database and schema
def go_db():
    trace(1, "Initialising database...")
    dbo = DB()
    dbo.create_schema()
    return dbo

# call out to bash script to parse the vcf file
# todo - pull this under python source?
def getVCFvariants(FILE):
    from subprocess import Popen, PIPE, STDOUT
    cmd = os.path.join(config['REDUX_ENV'], 'getVCFvariants.sh')

    p = Popen([cmd, '-'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    p_stdout = p.communicate(input=FILE.read())[0]
    output = p_stdout.decode('utf-8')
    # trace(1, 'command output: {}'.format(output))
    return output

#routines - "arghandler" (sort prototype) - Zak

# I'm not sure what this procedure is for
def go_sort_db():
    #trace(0,"** process SNP data.")
    dbo = go_db()
    dbo.insert_sample_sort_data()
    #dbo.commit()
    dbo.sort_data()
    #trace(0,"** + SNP processing done.")

    
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
        vcf_calls = analyzeVcf(config['unzip_dir'] + file)
        bed_calls = analyzeBed(config['unzip_dir'] + file)
        bed_index = [0]
        for lineno in range(len(d)):
            d[lineno] += ','
            if s[lineno] in vcf_calls:
                d[lineno] += vcf_calls[s[lineno]]
            d[lineno] += makeCall(s[lineno], bed_index, bed_calls)

        for line in d:
            print (line)

# routines - Iain 


# Returns a dict of position -> mutation mappings
# Modified from Harald's analyzeVCF, this version returns every mutation with
# its derived value, regardless of whether it was ancestral or not
def readHg19Vcf(file):

    with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
        trace (30, "   Extracting VCF: %s" % vcffile)
        result = {}
        for line in vcffile:
            fields = line.split()
            if (fields[0] == 'chrY' and int(fields[1]) > 0 and fields[3] != '.'
                    and fields[4] != '.'):
                result[fields[1]] = [int(fields[1]), str(fields[3]), str(fields[4])]
        return result


# get number of lines in a file
# probably will fail on utf-8 encoding or binary files in general
def file_len(fname):
    #File length, thanks to StackOverflow
    #https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    i=-1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# populate agebed table from age.bed
# This procedure dumps what may be in the agebed table and replaces it with the
# ranges defined in age.bed
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
        dbo.dc.execute('delete from agebed')
        dbo.dc.execute('drop table if exists tmpt')
        dbo.dc.execute('create temporary table tmpt(a,b,c)')
        dbo.dc.executemany('insert into tmpt values(?,?,?)', ranges)
        dbo.dc.execute('''insert or ignore into bedranges(minaddr,maxaddr)
                           select b,c from tmpt''')
        dbo.dc.execute('''insert into agebed(bID)
                           select id from bedranges
                           inner join tmpt t on t.b=minaddr and t.c=maxaddr''')
        dbo.dc.execute('drop table tmpt')


# populate a table of STR definitions
# Ordering is optional. If not given, table stores SNP names in FTDNA Y111 order.
# Otherwise, It's a 111-number vector.
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
    if ordering and (len(ordering)==len(strdefs)):
        for tup in zip(ordering,strdefs):
            dbo.dc.execute('insert or ignore into strs(ordering,strname) values(?,?)', tup)
    else:
        for tup in enumerate(strdefs):
            dbo.dc.execute('insert or ignore into strs(ordering,strname) values(?,?)', tup)


# pull information about the kits from the web api
# if API==None, read from
def get_kits (API='http://haplogroup-r.org/api/v1/uploads.php', qry='format=json'):
    import requests, json
    try:
        # choose where to pull the kit data
        if not API:
            trace(1, 'reading kit info from json.out')
            js = json.loads(open('json.out').read())
        else:
            trace(1, 'reading kit info from the web')
            url = '?'.join([API, qry])
            res = requests.get(url)
            print('Res.encoding is', res.encoding)
            js = res.json()
            open('json.out','w').write(json.dumps(js))
    except:
        print('Failed to pull kit metadata from {}'.format(API))
        raise # fixme - what to do on error
        
    return js

# update the information about the kits
# the input is js, a json object that comes from the Haplogroup-R DW API
def update_metadata(db, js):
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
    trace(1, 'first row out of {}:{}'.format(len(rows),rows[0]))
    trace(1, '{} unique kit ids'.format(len(set([v['kitId'] for v in js]))))
    trace(1, '{} null surname'.format(len([v['surname'] for v in js if not v['surname']])))

    # populate the dependency tables
    # (testtype,isNGS) goes into testtypes
    tups = [y for y in set([(r[-2],r[-1]) for r in rows])]
    db.dc.executemany('insert or ignore into testtype(testNm,isNGS) values(?,?)', tups)

    # (surname+kitId+build) goes into person
    # fixme: we need a DNA-to-person mapping. This is a big kludge
    tups = [y for y in set([(r[-3],r[0],r[-4]) for r in rows])]
    db.dc.executemany('insert or ignore into person(surname,firstname,middlename) values(?,?,?)',
                          tups)

    for tbl,val,idx in (('country','country',-7),
                        ('surname', 'surname',-3),
                        ('origin', 'origin',-6),
                        ('lab', 'labNm',-5),
                        ('build', 'buildNm',-4)):
        tups = [(y,) for y in set([v[idx] for v in rows])]
        trace(1, 'first tuple to insert into {}: {}'.format(tbl, tups[0]))
        db.dc.executemany('insert or ignore into {}({}) values(?)'.format(tbl,val),
                              tups)

    # create temporary table, where columns correspond to values above
    db.dc.execute('''create temporary table tmpt(
           a TEXT, b TEXT, c TEXT, d TEXT,
           e TEXT, f TEXT, g TEXT, h TEXT,
           i TEXT,
           j TEXT, k TEXT, l TEXT, m TEXT,
           n TEXT, o TEXT, p TEXT)''')
    db.dc.executemany('''INSERT INTO tmpt(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', rows)
    db.dc.execute('''
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
    trace(1, '{} rows inserted into dataset'.format(db.dc.execute('select count(*) from dataset').fetchone()[0]))
    db.dc.execute('drop table tmpt')

# load data into the contig table
def populate_contigs(db):
    trace(1,'populate contigs table')
    bid = get_build_byname(db, 'hg19')
    db.dc.execute('insert into Contig(buildID,description,length) values(?,?,?)',
                      (bid, 'chrY', 59373566))
    bid = get_build_byname(db, 'hg38')
    db.dc.execute('insert into Contig(buildID,description,length) values(?,?,?)',
                      (bid, 'chrY', 57227415))

# populate dataset information from Haplogroup-R data warehouse API
def populate_fileinfo(dbo, fromweb=True):
    if fromweb:
        js = get_kits()
    else:
        js = get_kits(API=None)
    update_metadata(dbo, js)

# update snp definitions from a csv DictReader instance
# fixme - update snpnames
def updatesnps(db, snp_reference, buildname='hg38'):
    bid = get_build_byname(db, buildname)
    db.dc.execute('drop table if exists tmpt')
    db.dc.execute('create temporary table tmpt(a integer, b integer, c text, d text, e text, unique(a,b,c,d,e))')
    db.dc.executemany('INSERT OR IGNORE INTO tmpt(a,b,c,d,e) VALUES (?,?,?,?,?)',
             ((bid, rec['start'],
                   rec['allele_anc'].strip(), rec['allele_der'].strip(), rec['Name'])
                     for rec in snp_reference))
    db.dc.execute('''insert or ignore into alleles(allele) 
                        select distinct c from tmpt''')
    db.dc.execute('''insert or ignore into alleles(allele) 
                        select distinct d from tmpt''')
    db.dc.execute('''insert or ignore into variants(buildID, pos, anc, der)
                        select a, b, an.id, dr.id from tmpt
                        inner join alleles an on an.allele = c
                        inner join alleles dr on dr.allele = d''')
    db.dc.execute('''insert or ignore into snpnames(snpname,vID)
                        select t.e, v.id from tmpt t, variants v, alleles a, alleles d
                        where t.c=a.allele and v.anc=a.id and
                              t.d=d.allele and v.der=d.id and
                              v.pos = t.b and v.buildID=t.a''')

    db.dc.execute('drop table tmpt')

# pull SNP definitions from the web at ybrowse.org
# this should be called after the build table is populated
# refresh files if they are older than maxage; do nothing if maxage < 0
def get_SNPdefs_fromweb(db, maxage, url='http://ybrowse.org/gbrowse2/gff'):
    import urllib, time
    if maxage < 0:
        return
    # convert to seconds
    maxage = 24*3600*maxage
    for (build,) in db.dc.execute('select buildNm from build'):
        deltat = maxage + 1
        fbase = 'snps_{}.csv'.format(build)
        fget = os.path.join(url, fbase)
        fname = os.path.join(config['REDUX_DATA'], fbase)
        try:
            if os.path.exists(fname):
                deltat = time.clock() - os.path.getmtime(fname)
            if deltat > maxage:
                trace (1, 'refresh: {}'.format(fbase))
                urllib.request.urlretrieve(fget, fname)
                deltat = time.clock() - os.path.getmtime(fname)
        except:
            pass
        if not os.path.exists(fname) or deltat > maxage:
            trace(0, 'failed to update {} from the web'.format(fname))
    return

# populate SNP definitions; refresh from web if we have is older than maxage
# (in days)
def populate_SNPs(dbo, maxage=config['max_snpdef_age']):
    get_SNPdefs_fromweb(dbo, maxage=maxage)
    # update known snps for hg19 and hg38
    with open(os.path.join(config['REDUX_DATA'], config['b37_snp_file'])) as snpfile:
        snp_reference = csv.DictReader(snpfile)
        updatesnps(dbo, snp_reference, 'hg19')
    with open(os.path.join(config['REDUX_DATA'], config['b38_snp_file'])) as snpfile:
        snp_reference = csv.DictReader(snpfile)
        updatesnps(dbo, snp_reference, 'hg38')

# return a) sum of BED ranges for kit (how many bases are covered by the test)
# and b) sum of coverage that is in the age BED range. These stats are needed
# for age calculations. This is not done in the "brute force" way because it
# can be compute intensive to search the list for every range.
def get_kit_coverage(dbo, pid):
    br = dbo.dc.execute('''select 1,minaddr from bedranges r
                           inner join bed b on b.bid=r.id and b.pid=?
                                union
                           select 1,maxaddr from bedranges r
                           inner join bed b on b.bid=r.id and b.pid=?
                                union
                           select 2,minaddr from bedranges b, agebed a
                           where a.bID=b.id
                                union
                           select 2,maxaddr from bedranges b, agebed a
                           where a.bID=b.id
                           order by 2,1''', (pid,pid))
    ids = {1:0, 2:1}
    accum1 = 0
    toggles = [False, False]
    post = None
    for r in br:
        toggles[ids[r[0]]] ^= True
        if post and (toggles[0] ^ toggles[1]):
            accum1 += r[1]-post
            post = None
        elif toggles[0] & toggles[1]:
            post = r[1]
    if post:
        accum1 += r[1]-post
    accum2 = dbo.dc.execute('''select sum(maxaddr-minaddr) from bedranges r
                           inner join bed b on r.id=b.bid and b.pid=?''',
                        (pid,)).fetchone()[0]
    return accum2, accum1


# efficiently determine if a set of positions is contained in a set of ranges
# return vector of True values if v contained in a range, else False
# does not consider the endpoints
# ranges must be sorted on their minaddr
# input vector must be sorted
def in_range(v_vect, ranges, spans):
    c_vect = []
    ii = 0
    nmax = len(ranges)
    for iv,v in enumerate(v_vect):
        while ii < nmax and ranges[ii][0] < v:
            ii += 1
        if ii > 0 and v > ranges[ii-1][0] and v < ranges[ii-1][1]:
            if spans and v+spans[iv] < ranges[ii-1][1]:
                # span fits within the range
                c_vect.append(True)
            elif spans:
                # span exceeded the upper end of the range
                c_vect.append(False)
            else:
                # no span - SNP fits within range
                c_vect.append(True)
        else:
            c_vect.append(False)
    if len(c_vect) != len(v_vect):
        raise ValueError
    return c_vect

# Check kit coverage for a vector of variants; return a coverage vector that
# indicates if that person has a BED range such that the call is a) in a range,
# b) on the lower edge of a range, c) on the upper edge of a range, or d) not
# covered by a range. If spans is passed, it corresponds to the maximum length
# affected by the variant in vids: 1 for SNPs, max(ref,alt) for indels
def get_call_coverage(dbo, pid, vids, spans=None):
    # FIXME currently only returns True/False, doesn't handle range ends
    # FIXME maybe more efficient to pass positions instead of ids
    dc = dbo.cursor()
    dc.execute('drop table if exists tmpt')
    dc.execute('create table tmpt (vid integer)')
    trace(3, 'get_call_coverage: vids:{}...'.format(vids[:10]))
    if spans:
        trace(3, 'get_call_coverage: spans:{}...'.format(spans[:10]))
    dc.executemany('insert into tmpt values(?)', [(a,) for a in vids])
    calls = dc.execute('''select v.id, v.pos from variants v
                          inner join tmpt t on t.vid=v.id
                          order by 2''')
    # form a list of positions of interest
    cv = [v[1] for v in calls]
    trace(500, '{} calls: {}...'.format(len(cv), cv[:20]))
    rc = dbo.cursor()
    ranges = rc.execute('''select minaddr,maxaddr from bedranges r
                           inner join bed b on b.bID=r.id
                           where b.pid=?
                           order by 1''', (pid,))
    # form a list of ranges of interest
    rv = [v for v in ranges]
    if len(rv) == 0:
        return []
    trace(500, '{} ranges: {}...'.format(len(rv), rv[:20]))
    coverage = in_range(cv, rv, spans)
    dc.close()
    rc.close()
    return coverage

# populate regions from a FTDNA BED file
# fname is an unpacked BED file
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

# experimental - pack call information into an integer
# pos, anc, der, passfail, q1, q2, nreads, passrate
def pack_call(call_tup):
    maxshift = 31
    nqbits = 7
    ncbits = 9
    npbits = 8
    # pack in a pass/fail flag
    bitfield = {'PASS': 1<<maxshift, 'FAIL': 0}[call_tup[3]]
    maxshift -= 1
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

# experimental - unpack that corresponds to pack_call
# if pack_call changes, this needs to change
def unpack_call(bitfield):
    maxshift = 31
    nqbits = 7
    ncbits = 9
    npbits = 8
    passrate = float(bitfield & (1<<npbits) - 1) / ((1<<npbits) - 1)
    calls = (bitfield>>npbits) & ((1<<ncbits) - 1)
    q2 = float((bitfield>>(npbits+ncbits)) & ((1<<nqbits) - 1)) / 2.
    q1 = float((bitfield>>(npbits+ncbits+nqbits)) & ((1<<nqbits) - 1)) /3.
    if (1<<maxshift) & bitfield:
        passfail = True
    else:
        passfail = False
    return passfail, q1, q2, calls, passrate

# populate calls, quality, and variants from a VCF file
# fname is an unzipped VCF file
def populate_from_VCF_file(dbo, bid, pid, fileobj):
    dc = dbo.cursor()
    b = dc.execute('select buildNm from build where id=?', (bid,)).fetchone()[0]
    if b != 'hg38':
        trace(0, 'currently unable to parse build {}'.format(b))
        return
    # parse output of getVCFvariants(fname)
    parsed = getVCFvariants(fileobj)
    tups = []
    for line in parsed.splitlines():
        # pos, anc, der, passfail, q1, q2, nreads, passrate
        tups.append(line.split())

    # trace(1, 'parsed vcf: {}...'.format(tups[:3]))
    # filter down to the calls we want to store
    # fixme this is probably not the correct filtering
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
    dc.execute('''create temporary table tmpt(a integer, b integer,
                  c text, d text, e integer, f integer)''')
    dc.executemany('insert into tmpt values(?,?,?,?,?,?)',
                          [[bid]+v[0:3]+[pid]+[v[-1]] for v in call_info])
    # fixme - performance
    trace(3,'VCF update variants at {}'.format(time.clock()))
    dc.execute('''insert or ignore into variants(buildID, pos, anc, der)
                  select a, b, an.id, dr.id from tmpt
                  inner join alleles an on an.allele = c
                  inner join alleles dr on dr.allele = d''')
    trace(3,'done at {}'.format(time.clock()))
    trace(3,'VCF update calls at {}'.format(time.clock()))
    dc.execute('''insert into vcfcalls (pid,vid,callinfo)
                  select e, v.id, f from tmpt
                  inner join alleles an on an.allele = c
                  inner join alleles dr on dr.allele = d
                  inner join variants v on v.buildID = a and v.pos = b
                  and v.anc=an.id and v.der=dr.id''')
    trace(3,'done at {}'.format(time.clock()))
    dc.execute('drop table tmpt')
    dc.close()

    trace(500, '{} vcf calls: {}...'.format(len(passes), passes[:5]))
    return

# unpack any zip from FTDNA that has the bed file and vcf file
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

# unpack all zip files that can be handled from the HaplogroupR catalog If the
# zip file exists, unpack it. Future - skip if already loaded?  Walk through
# dataset table to find these files. Assume we've already pulled the list from
# the H-R web API and downloaded zip files. Future - download file?
def populate_from_dataset(dbo):
    import zipfile, re
    trace(1, 'populate from dataset with kit limit {}'.format(config['kitlimit']))
    dc = dbo.cursor()
    bed_re = re.compile(r'(\b(?:\w*[^_/])?regions(?:\[\d\])?\.bed)')
    vcf_re = re.compile(r'(\b(?:\w*[^_/])?variants(?:\[\d\])?\.vcf)')
    dc = dc.execute('select fileNm,buildID,ID from dataset')
    allsets = list([(t[0],t[1],t[2]) for t in dc])
    pc = dbo.cursor()
    pl = pc.execute('select distinct pid from vcfcalls')
    pexists = [p[0] for p in pl]
    trace(5,'allsets: {}'.format(allsets[:config['kitlimit']]))
    nkits = 0

    # fixme - hack - better index handling (drop index for insert performance)
    trace(3, 'drop indexes at {}'.format(time.clock()))
    dc.execute('drop index bedidx')
    dc.execute('drop index vcfidx')
    dc.execute('drop index vcfpidx')
    trace(3, 'done at {}'.format(time.clock()))

    for (fn,buildid,pid) in allsets:
        # if there are already calls for this person, skip the load
        # fixme - should also skip if BED entries exist
        if (not config['drop_tables']) and pid in pexists:
            trace(2, 'calls exist - skip {}'.format(fn[:50]))
            continue
        zipf = os.path.join(data_path('HaplogroupR'), fn)
        if not os.path.exists(zipf):
            trace(10, 'not present: {}'.format(zipf))
            continue
        try:
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
                    trace(2, 'populate from bed {}'.format(bedfile))
                    populate_from_BED_file(dbo, pid, bedf)
                with zf.open(vcffile,'r') as vcff:
                    trace(2, 'populate from vcf {}'.format(vcffile))
                    populate_from_VCF_file(dbo, buildid, pid, vcff)
            nkits += 1
        except:
            trace(0, 'FAIL on file {} (not loaded)'.format(zipf))
            # raise
        if nkits >= config['kitlimit']:
            break
        # fixme - if BED passes and VCF fails, cruft is left behind. This loop
        # could go into a transaction, but that slows the loading. I can't
        # think of any harm storing the extra BED info.

    # fixme - hack
    trace(3, 're-create indexes at {}'.format(time.clock()))
    dc.execute('create index bedidx on bed(pID,bID)')
    dc.execute('create index vcfidx on vcfcalls(vID)')
    dc.execute('create index vcfpidx on vcfcalls(pID)')
    dc.close()
    trace(3, 'done at {}'.format(time.clock()))

    # return without calculating coverage for large number of kits
    if config['kitlimit'] > 20:
        return 

    trace(1, 'calculate coverages')
    # below: maybe not part of populating dataset?
    # get variants in ascending order of number of calls across kits
    vc = dbo.cursor().execute('''select v.id,v.anc,v.der,count(c.pid)
                           from variants v
                           inner join vcfcalls c on c.vid = v.id
                           group by 1,2,3 order by 4''')
    # variant list is used more than once, so make a list
    vl = list([s[0] for s in vc if s[3] > 1])
    vc.close()
    trace(2,'closed cursor')
    # fixme - need to do similar for indels (save/reuse the list)
    fl = dbo.cursor()
    fl = fl.execute('''select distinct fileNm,buildID,ID from dataset d
                       inner join vcfcalls c on c.pid=d.id''')
    allsets = list([(t[0],t[1],t[2]) for t in fl])
    fl.close()
    try:
        for (fn,buildid,pid) in allsets:
            trace(1, 'calculate coverage for id {}'.format(pid))
            coverage = get_call_coverage(dbo, pid, vl)
            trace(1, 'calculate indel coverage for id {}'.format(pid))
            indel_coverage = get_indel_coverage(dbo, pid)
            if len(coverage) > 0:
                trace(2, 'coverage vector: {}...'.format(coverage[:20]))
            else:
                trace(2, 'no coverage calculated for {}'.format(pid))
    except:
        trace(0, 'did not calculate coverage')
        raise
    return


# initial database creation and table loads
def db_creation():
    db = DB(drop=config['drop_tables'])
    if config['drop_tables']:
        db.create_schema()
    populate_fileinfo(db, fromweb=config['use_web_api'])
    populate_STRs(db)
    populate_SNPs(db)
    populate_contigs(db)
    populate_age(db)
    return db
