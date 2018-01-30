#!/usr/bin/env python3
# coding: utf-8

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand, Zak Jones
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# run this script as a command with no args to execute some tests


import os,yaml,shutil,glob,re,csv,zipfile,subprocess
from db import DB
from collections import defaultdict
from array_api import *


REDUX_CONF = 'config.yaml'
config = yaml.load(open(REDUX_CONF))


# add src and bin directories to path
import sys
sys.path.insert(0, config['REDUX_ENV'])
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

# I don't know what this procedure is for; it looks like it needs to go in the
# main loop - I've removed some code for vcf files, which is handled
# differently
def skip_to_Hg19(dbo):

    # skip to <= 1 - unpack zips

    # skip to <= 10 - associate kits with people

    # skip to <= 11 - generate dictionary of variant positions 

    #   The VCF files are parsed twice. The first pass identifies the list
    #   of variants to be queried. The second pass reads the calls for those
    #   variants. This means we only need to treat positions with a variant,
    #   rather than every position in a chromosome.
    #   Using a dictionary here means only one copy of each variant is
    #   created.

    if (config['skip_to'] <= 11):

        #vcffiles

        trace (2, "Generating database of all variants...")
        
        #variant_dict

        #print(REDUX_ENV)
        variant_dict = {}

        # dump variant dict into sorted array

        trace (20, "   Dumping variants into array...")
        variant_array = np.array(list(variant_dict.values()))

        # variant_array = np.array([],dtype={'names': ('start', 'anc', 'der'),'formats': ('i4', 'S20', 'S20')})

        trace (30, "      Check variant [0] is %s" % variant_array[0])
        trace (30, "      Check variant [0] position is %s" % variant_array[0][1])
        trace (30, "      Check variant [%s] is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1]))
        trace (30, "      Check variant [%s] position is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1][1]))

        #db calls

        trace (20, "   Inserting data into variant array database...")
        dbo.insert_variants(variant_array, 'hg38')
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
        
    # skip to <= 12 - reading calls for variants
    
    #db calls
    
    # skip to <= 13 - name variants and derive ancestral values

    # Some variants are positive in the reference sequence, so we need to
    # look up their ancestral values. We'll get the SNP names while we're
    # at it.

    #db calls

    if (config['skip_to'] <= 13):
        # Read in SNPs from reference lists
        trace (2, "Getting names of variants...")
        trace (10, "   Importing SNP reference lists...")


        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
            
    # Print final message and exit {{{{

    t = float((time.clock() - start_time))
    trace (1, "Execution finished in: %.3f seconds" % t)


# routines - arghandler (redux1) - Zak

def go_all():
    go_backup()
    go_prep()
    go_db()
    
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
    

# I'm not sure what this procedure is for ???
def go_prep():

    trace(0,"** prepare file structure.")

    # SKIPZIP check (beg)

    if not config['skip_zip']:

        # Check ZIPDIR - contains existing zip files

        if not os.path.exists(config['zip_dir']):
            trace(0,"Input zip folder does not appear to exist. Aborting.\n")
            sys.exit()

        # Check WORKING - the zip working exists && Empty it, otherwise make it

        refresh_dir('working')

        # Check UNZIPDIR - contains zip output; refresh it

        refresh_dir('unzips',not config['zip_update_only'])

        # Get the list of input files

        FILES = glob.glob(os.path.join(config['REDUX_ENV'], 'zips', 'bigy-*.zip'))

        if len(FILES) == 0:
            trace(0,"No input files detected in zip folder. Aborting.")
            trace(0,"Check name format: should be bigy-<NAME>-<NUMBER>.zip\n")
            sys.exit()
        else:
            trace(0,"input files detected: " + str(FILES))

        # Check whether unzip is installed

        if not cmd_exists('unzip'):
            trace(0,"Unzip package not found. Aborting.")
            sys.exit()

        # Check whether SNP list exists - handled elsewhere

        # Check whether merge-ignore list exists

        touch_file('merge-ignore.txt')

        # fix bash code (beg)

        # Unzip each zip in turn - handled elsewhere

    # fix bash code (end)

    trace(0,"** + prep done.")
    

# routines - arghandler (redux2) - Zak

# utility procuedure to drop/create the database and schema
def go_db():
    trace(1, "Initialising database...")
    dbo = DB()
    dbo.create_schema()
    return dbo

#note: it looks like clades.py is doing something like this for hg19
def getH38references():
    foo = 1

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

#routines - "arghandler" (new v2 schema)- Jef/Zak

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
            dbo.dc.execute('insert into strs(ordering,strname) values(?,?)', tup)
    else:
        for tup in enumerate(strdefs):
            dbo.dc.execute('insert into strs(ordering,strname) values(?,?)', tup)


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
            js = json.loads(res.content)
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
    db.dc.executemany('insert into testtype(testNm,isNGS) values(?,?)', tups)

    # (surname+kitId+build) goes into person
    # fixme: we need a DNA-to-person mapping. This is a big kludge
    tups = [y for y in set([(r[-3],r[0],r[-4]) for r in rows])]
    db.dc.executemany('insert into person(surname,firstname,middlename) values(?,?,?)',
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
        INSERT INTO dataset(kitId, importDt, fileNm, lng,
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
    for (build,) in db.dc.execute('select buildNm from build'):
        deltat = maxage + 1
        fbase = 'snps_{}.csv'.format(build)
        trace (1, 'refresh: {}'.format(fbase))
        fget = os.path.join(url, fbase)
        fname = os.path.join(config['REDUX_DATA'], fbase)
        try:
            if os.path.exists(fname):
                deltat = time.time() - os.path.getmtime(fname)
            if deltat > maxage:
                urllib.request.urlretrieve(fget, fname)
                deltat = time.time() - os.path.getmtime(fname)
        except:
            pass
        if not os.path.exists(fname) or deltat > maxage:
            trace(0, 'failed to update {} from the web'.format(fname))
    return

# populate SNP definitions; refresh from web if we have is older than maxage (seconds)
def populate_SNPs(dbo, maxage=3600*24*5):
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
def in_range(v_vect, ranges):
    c_vect = []
    ii = 0
    nmax = len(ranges)
    for v in v_vect:
        while ii < nmax and ranges[ii][0] < v:
            ii += 1
        if ii > 0 and v > ranges[ii-1][0] and v < ranges[ii-1][1]:
            c_vect.append(True)
        else:
            c_vect.append(False)
    if len(c_vect) != len(v_vect):
        raise ValueError
    return c_vect


# check kit coverage for a vector of variants return a vector of BED range
# coverage for a person that represents if the call is a) in a range, b) on the
# lower edge of a range, c) on the upper edge of a range, or d) not covered by
# a range
def get_call_coverage(dbo, pid, vids):
    # FIXME currently only returns True/False, doesn't handle range ends
    dc = dbo.cursor()
    dc.execute('drop table if exists tmpt')
    dc.execute('create table tmpt (vid integer)')
    dc.executemany('insert into tmpt values(?)', [(a,) for a in vids])
    calls = dc.execute('''select v.id, v.pos from variants v
                          inner join tmpt t on t.vid=v.id
                          order by 2''')
    cv = [v[1] for v in calls]
    trace(500, '{} calls: {}...'.format(len(cv), cv[:20]))
    rc = dbo.cursor()
    ranges = rc.execute('''select minaddr,maxaddr from bedranges r
                           inner join bed b on b.bID=r.id
                           where b.pid=?
                           order by 1''', (pid,))
    rv = [v for v in ranges]
    if len(rv) == 0:
        return []
    trace(500, '{} ranges: {}...'.format(len(rv), rv[:20]))
    coverage = in_range(cv, rv)
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
        passes = [t for t in tups if t[2] != '.' and (t[3] == 'PASS' or
                  (int(t[6]) < 4 and float(t[7]) > .75))]
    except ValueError:
        trace(0, 'parsing VCF failed on {}'.format(t))

    # save the distinct alleles
    alleles = set([x[1] for x in passes] + [x[2] for x in passes])
    dc.executemany('insert or ignore into alleles(allele) values(?)',
                       [(x,) for x in alleles])

    # save the call quality info
    # fixme - update schema to handle quality info
    pass

    # execute sql on results to save in vcfcalls
    dc.execute('drop table if exists tmpt')
    dc.execute('''create temporary table tmpt(a integer, b integer,
                  c text, d text, e integer)''')
    dc.executemany('insert into tmpt values(?,?,?,?,?)',
                          [[bid]+v[0:3]+[pid] for v in passes])
    dc.execute('''insert or ignore into variants(buildID, pos, anc, der)
                  select a, b, an.id, dr.id from tmpt
                  inner join alleles an on an.allele = c
                  inner join alleles dr on dr.allele = d''')
    dc.execute('''insert into vcfcalls (pid,vid)
                  select e, v.id from tmpt
                  inner join alleles an on an.allele = c
                  inner join alleles dr on dr.allele = d
                  inner join variants v on v.buildID = a and v.pos = b
                  and v.anc=an.id and v.der=dr.id''')
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
    fl = dc.execute('select fileNm,buildID,DNAID from dataset')
    allsets = list([(t[0],t[1],t[2]) for t in fl])
    trace(1,'allsets: {}'.format(allsets[:config['kitlimit']]))
    nkits = 0
    for (fn,buildid,dnaid) in allsets:
        zipf = os.path.join(data_path('HaplogroupR'), fn)
        if not os.path.exists(zipf):
            # trace(2, 'not present: {}'.format(zipf))
            continue
        with zipfile.ZipFile(zipf) as zf:
            trace(1, 'populate from {}'.format(zf))
            listfiles = zf.namelist()
            bedfile = vcffile = None
            for ff in listfiles:
                dirname, basename = os.path.split(ff)
                if bed_re.search(basename):
                    bedfile = ff
                elif vcf_re.search(basename):
                    vcffile = ff
            if (not bedfile) or (not vcffile):
                trace(0, 'WARN: missing data in '+zipf)
                continue
            with zf.open(bedfile,'r') as bedf:
                trace(1, 'populate from bed {}'.format(bedfile))
                populate_from_BED_file(dbo, dnaid, bedf)
            with zf.open(vcffile,'r') as vcff:
                trace(1, 'populate from vcf {}'.format(vcffile))
                populate_from_VCF_file(dbo, buildid, dnaid, vcff)
        nkits += 1
        if nkits > config['kitlimit']:
            break

    dc.close()
    # below: maybe not part of populating dataset?
    # get variants in ascending order of number of calls across kits
    vc = dbo.cursor().execute('''select v.id,v.anc,v.der,count(c.pid)
                           from variants v
                           inner join vcfcalls c on c.vid = v.id
                           group by 1,2,3 order by 4''')
    # variant list if called more than once
    vl = list([s[0] for s in vc if s[3] > 1])
    for (fn,buildid,dnaid) in allsets:
        coverage = get_call_coverage(dbo, dnaid, vl)
        if len(coverage) > 0:
            trace(1, 'coverage vector: {}...'.format(coverage[:20]))


# initial database creation and table loads
def db_creation():
    db = DB(drop=True)
    db.create_schema()
    populate_fileinfo(db, fromweb=False)
    populate_STRs(db)
    populate_SNPs(db)
    populate_contigs(db)
    populate_age(db)
    return db


# test framework
# The following code is not really part of the program; it's for calling
# procedures while developing or testing.
if __name__ == '__main__':
    # print (config)
    db = db_creation()
    populate_from_dataset(db)
    ids = get_dna_ids(db)
    out = get_variant_csv(db,ids)
    open('csv.out','w').write(out)
    db.commit()
    db.close()

sanity_tests='''
./lib.py
should produce no errors and produce variants.db after some time

sqlite3 variants.db
sqlite> .mode tabs
sqlite> .header on

sqlite> select vid, snpname, der, allele from snpnames inner join variants v on v.id=vid and v.pos=6753258 inner join alleles a on v.der=a.id;

should return nine rows similar to
vID	snpname	der	allele
7574	L147.2	3	C
7574	L147.5	3	C
7574	L147.3	3	C
179549	L1283	4	G
181076	Y8884	2	A
7574	L147.4	3	C
7574	PF4883	3	C
7574	L147.1	3	C
7574	L147	3	C

sqlite> select count(*) from variants;
613049

see some other sample queries in ../examples directory
'''
