#!/usr/bin/env python3
"""
  Purpose:
    proof of concept stats and other operations on bed files from BigY tests
    since it's only a demo, nothing has been optimized, and calculation may not
    be quite correct

  Usage:
    edit file & run; expects an "unzip" directory
    clades.py create : populate database from "unzip"; prereq for other ops
    clades.py

  Copyright:
    For free distribution under the terms of the
    GNU General Public License, version 3 (29 June 2007)
    https://www.gnu.org/licenses/gpl.html

  Jef Treece, 22 Mar 2017
"""

import sqlite3, os, time, sys
from collections import defaultdict
import random
import argparse
import locale
import csv
import subprocess

# affects diagnostic messages, which always go to stderr
# this is only the default; increase with -v flag(s)
DEBUG=0

# where the files reside, currently hardwired
UNZIPDIR = 'unzip'

# the name of the SQLite database on disk
DBFILE = 'clades.db'

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='count')
parser.add_argument('action', nargs='*')
parser.add_argument('-s', '--snp', nargs=1)
parser.add_argument('-i', '--implications', action='store_true')
parser.add_argument('-t', '--tree', action='store_true')
parser.add_argument('-b', '--badlist', action='store_true')
parser.add_argument('-k', '--kits', action='store_true')
namespace = parser.parse_args(sys.argv[1:])
verbose = vars(namespace)['verbose']
if not verbose:
    verbose = DEBUG

if namespace.snp:
    querysnp = vars(namespace)['snp'][0]
else:
    querysnp = False

# create: new database from all of the .bed and vcf files
# fastest if you remove the old .db file first
create = stats1 = stats2 = docalls = listfiles = listbed = mergeup = updatesnps = False
for a in namespace.action:
    if a == 'create':
        create = True
    elif a == 'stats1':
        stats1 = True
    elif a == 'stats2':
        stats2 = True
    elif a == 'docalls':
        docalls = True
    elif a == 'listfiles':
        listfiles = True
    elif a == 'listbed':
        listbed = True
    elif a == 'updatesnps':
        updatesnps = True
    elif a == 'mergeup':
        mergeup = True
    else:
        print('unknown:', action, 'exiting')

# trace (print) output if it exceeds the noise level
def trace (level, msg):
    if level <= verbose:
        print(msg, file=sys.stderr)
        sys.stderr.flush()

tabledefs = [
    # for "pseudo" files that did not originate from test data, kit is 0
    'files(ID INTEGER PRIMARY KEY, name TEXT, seq INTEGER DEFAULT 0, kit INTEGER DEFAULT 0)',
    'bed(ID INTEGER REFERENCES files(ID), minaddr INTEGER, maxaddr INTEGER)',
    # vcfcalls mostly comes from the .vcf files, so origin is "k"
    # there are some tweaks from implications.txt, and origin is set to "^" or "<"
    'vcfcalls(ID INTEGER REFERENCES files(ID), vid REFERENCES VARIANTS(id), origin TEXT DEFAULT "k")',
    'vcfstats(ID INTEGER REFERENCES files(ID), ny INTEGER, nv INTEGER, ns INTEGER, ni INTEGER, nr INTEGER)',
    'bedstats(ID INTEGER REFERENCES files(ID), coverage1 INTEGER, coverage2 INTEGER, nranges INTEGER)',
    'vcfrej(ID INTEGER REFERENCES files(ID), vid INTEGER REFERENCES variants(id), origin TEXT DEFAULT "k")',
    'meta(key TEXT, val TEXT)',
    'variants(ID INTEGER PRIMARY KEY, pos INTEGER, ref TEXT, alt TEXT, UNIQUE(pos, ref, alt))',
    'snpnames(ID INTEGER REFERENCES variants(ID), name TEXT)',
	]
indexdefs = [
    'rangeidx on bed(id, minaddr, maxaddr)',
    'rangeidx2 on bed(id, minaddr)',
    'fileidx on files(name)',
    'vcfidx on vcfcalls(id, vid)',
    'vcfidx2 on vcfcalls(id)',
    'vcfridx on vcfrej(id,vid)',
    'varidx on variants(pos)',
	]

t0 = time.time()

if create:
    try:
        os.unlink(DBFILE)
    except:
        pass

# create a sqlite database, cursor, and tables
dbconn = sqlite3.connect(DBFILE)
dbcurs = dbconn.cursor()

if create:
    updatesnps = True

    # nranges, coverage1 (total r1), coverage2 (gated by r2)
    def get_coverage(r1, r2):
        ids = {r1[0][0]:0, r2[0][0]:1}
        accum1 = 0
        toggles = [False, False]
        intervals = [(a,b) for (a,b,c) in r1+r2] + [(a,c) for (a,b,c) in r1+r2]
        post = None
        for r in sorted(intervals, key=lambda a: a[1]):
            toggles[ids[r[0]]] ^= True
            if post and (toggles[0] ^ toggles[1]):
                accum1 += r[1]-post
                post = None
            elif toggles[0] & toggles[1]:
                post = r[1]
        if post:
            accum1 += r[1]-post
        accum2 = 0
        for id, mina, maxa in r1:
            accum2 += maxa-mina
        return accum2, accum1

    for table in tabledefs:
        dbcurs.execute('CREATE TABLE %s' % table)
    for idx in indexdefs:
        dbcurs.execute('CREATE INDEX %s' % idx)

    # match order to the default output of ls
    # ls depends on such things as locale setting, aliases, etc
    # if in doubt, use the 'listfiles' action to produce a file list
    # this order matters if you need to match to columns of a spreadsheet
    locale.setlocale(locale.LC_ALL, '')
    files = sorted(os.listdir(UNZIPDIR), key=locale.strxfrm)

    # any metadata about the database we might want to know later
    meta = [('datadir', UNZIPDIR),
                ('date', ' '.join((time.asctime(), time.tzname[time.daylight]))),
                ('locale', 'LC_ALL'),
                ('directory', os.getcwd()),
                ]
    dbcurs.executemany('insert into meta values(?,?)', meta)


    ranges = []
    dbcurs.execute('INSERT INTO files(name) VALUES("age.bed")')
    myid = dbcurs.lastrowid
    lines = open('age.bed').readlines()
    for line in lines:
        chromo, minr, maxr = line.split()
        ranges.append((myid, int(minr), int(maxr)))
    dbcurs.executemany('insert into bed values(?,?,?)', ranges)
    cover_ranges = ranges

    trace(1, 'enter at {:.2f} seconds'.format(time.time() - t0))

    stat_list = []
    for nf,fname in enumerate(files):
        if fname.endswith('.bed'):
            dbcurs.execute('INSERT INTO files(name,kit,seq) VALUES(?,1,?)', (os.path.splitext(fname)[0],nf))
            myid = dbcurs.lastrowid
            ranges = []
            try:
                for line in open(os.path.join(UNZIPDIR, fname)):
                    ychr, minr, maxr = line.split()
                    ranges.append((myid, int(minr), int(maxr)))
            except:
                trace(0, 'FAILED on file {}'.format(fname))
                continue
            dbcurs.executemany('insert into bed values(?,?,?)', ranges)
            tcover, pcover = get_coverage(ranges, cover_ranges)
            stat_list.append((myid, tcover, pcover, len(ranges)))
    trace(1, 'bed done at {:.2f} seconds'.format(time.time() - t0))
    dbcurs.executemany('insert into bedstats values(?,?,?,?)', stat_list)


    # insert all of the special variants from implications.txt
    dbcurs.execute('INSERT INTO files(name) VALUES("implications")')
    myid = dbcurs.lastrowid
    ll = {}
    vl = []
    with open('implications.txt') as implications:
        for line in implications:
            toks = line.split()
            if toks and toks[0] in ('^', '<'):
                addr = toks[1]
                ref = toks[3]
                alt = toks[4]
                vl.append((addr, ref, alt))
                ll[(addr, ref, alt)] = toks[0]
    dbcurs.executemany('insert or ignore into variants(pos,ref,alt) values (?,?,?)', ll)
    cl = []
    for tup in vl:
        tid = dbcurs.execute('select id from variants where pos=? and ref=? and alt=?', tup).fetchone()[0]
        cl.append((myid, tid, ll[tup]))
    dbcurs.execute ('delete from vcfcalls where id=?', (myid,))
    trace(3, 'inserting to vcfcalls {}'.format(cl))
    dbcurs.executemany('insert into vcfcalls(id, vid, origin) values(?,?,?)', cl)

    calls = []
    rejects = []
    stats = []
    ipos = [a for (a,b,c) in vl]
    for fname in files:
        if fname.endswith('.vcf'):
            ny = nv = ns = ni = nr = 0
            myid = dbcurs.execute('select id from files where name = ?', (os.path.splitext(fname)[0],)).fetchone()[0]
            trace(3, 'reading file {} from {}'.format(fname, UNZIPDIR))
            with open(os.path.join(UNZIPDIR, fname)) as lines:
                for line in lines:
                    try:
                        chrom,pos,iid,ref,alt,qual,filt,inf,fmt,ii = line.split()
                    except:
                        continue
                    if chrom == 'chrY':
                        ny += 1
                        # keep certain rejects for swapping reference genome positions later
                        # rather than swap&call rejects, mark them dirty?
                        if filt == 'REJECTED':
                            nr += 1
                            if pos in ipos:
                                rejects.append((myid, pos, ref, alt))
                        elif filt == 'PASS':
                            nv += 1
                            if ref != '.' and alt != '.':
                                calls.append((myid, pos, ref, alt))
                                if len(ref) == 1 and len(alt) == 1:
                                    ns += 1 # snps
                                else:
                                    ni += 1 # indels
            trace(3, 'length of rejects now: {}'.format(len(rejects)))
            stats.append((myid, ny, nv, ns, ni, nr))
    dbcurs.executemany('INSERT INTO vcfstats(id,ny,nv,ns,ni,nr) VALUES(?,?,?,?,?,?)', stats)
    vl = []
    for iid, pos, ref, alt in rejects:
        try:
            vid = dbcurs.execute('select id from variants where pos=? and ref=? and alt=?', (pos,ref,alt)).fetchone()[0]
        except:
            dbcurs.execute('insert into variants(pos,ref,alt) values (?,?,?)', (pos,ref,alt))
            vid = dbcurs.lastrowid
        vl.append((iid,vid))
    dbcurs.executemany('INSERT INTO vcfrej(id,vid) VALUES(?,?)', vl)
    vl = []
    for iid, pos, ref, alt in calls:
        try:
            vid = dbcurs.execute('select id from variants where pos=? and ref=? and alt=?', (pos,ref,alt)).fetchone()[0]
        except:
            dbcurs.execute('insert into variants(pos,ref,alt) values (?,?,?)', (pos,ref,alt))
            vid = dbcurs.lastrowid
        vl.append((iid,vid))
    dbcurs.executemany('INSERT INTO vcfcalls(id,vid) VALUES(?,?)', vl)

    rejects = calls = stats = vl = None

    trace(1, 'vcf done at {:.2f} seconds'.format(time.time() - t0))
    dbconn.commit()
    trace(1, 'commit done at {:.2f} seconds'.format(time.time() - t0))



#####
# check for positive calls, covered, and not-covered in all kits
#
if docalls:

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

    c1 = dbconn.cursor()
    c2 = dbconn.cursor()

    # wishlist: redo the refpos swapping from implications
    # caution if adding new variants, might not capture all of the rejects
    # see above, in creating the database
    # redoing here permits quick edit of impl.txt & re-run without a re-build
    # best to do a re-build on final run

    # the reference assembly genome, hg19 was mainly someone U152+, so the
    # REF/ALT need to be swapped for some POSs these were read from the
    # implications file, a result of hand-correcting oddities this means people
    # who have a call for this POS in their vcf don't really have a SNP and
    # vice-versa these values are both in variants, below, and output is
    # corrected when we write at the end

    ref_swaps = {}
    swaps = c1.execute('select v.id,v.pos,v.ref,v.alt from variants v, vcfcalls c where c.vid=v.id and c.origin="<"')
    for vid,pos,ref,alt in swaps:
        # ref_swaps[vid] = (ref,alt) # as specified in implications
        trace(2,'swappos id={}:{}.{}.{}'.format(vid,pos,ref,alt))
        vc =c2.execute('select id from variants where pos=? and ref=? and alt=?', (pos,alt,ref)) # swapped version versus what's in implications
        for (vid,) in vc: # swapped version appears in variants
            trace(2,'swappos id={}:{}.{}.{}'.format(vid,pos,alt,ref))
            ref_swaps[vid] = (ref,alt)

    inserted = {}
    ins = c1.execute('select c.vid from vcfcalls c where c.origin="^"')
    for (vid,) in ins:
        trace(2,'inserted:variant id {}'.format(vid))
        inserted[vid] = None

    # list and count variants from the kits and uncount the inserted ones
    vcounts = {}
    vl = c1.execute('select count(id), vid from vcfcalls group by vid')
    for count, vid in vl:
        trace(3,'{}: {}'.format(vid,count))
        vcounts[vid] = count
    vl = c1.execute('select vid from vcfcalls where origin in ("^","<")')
    for (vid,) in vl:
        vcounts[vid] -= 1

    # which variants are SNPS
    snps = {}
    vl = c1.execute('select id from variants where length(ref)=1 and length(alt)=1')
    for (vid,) in vl:
        snps[vid] = True

    # which variants are indels
    indels = {}
    vl = c1.execute('select id from variants where (length(ref)>1 or length(alt)>1)')
    for (vid,) in vl:
        indels[vid] = True

    trace(2, '{} snps and {} indels'.format(len(snps), len(indels)))

    trace(1, 'variants gathered at {:.2f} seconds'.format(time.time() - t0))

    # algo
    # get all of the cbu, cbl, and cblu - this is fast
    # for each tester id
    #   get the calls for that id
    #   get the bed ranges for that id
    #   check coverage
    #      ';cbu' ';cbu' ';cblu'
    #      if not covered append ';nc'

    calls = defaultdict(dict)

    c1.execute('create temporary table tpos (pos int)')
    c1.execute('insert into tpos select distinct pos from variants')

    # these run fast - gather cbu,cbl,cblu all at once, and they do not depend
    # on REF,POS
    cbu = 'select id,maxaddr from bed where maxaddr in (select * from tpos)'
    cbl = 'select id,minaddr from bed where minaddr in (select * from tpos)'
    cblu = 'select id,maxaddr from bed where maxaddr=minaddr and maxaddr in (select * from tpos)'
    call_strings = ['', ';cblu', ';cbu', ';cbl', ';nc']

    for ii,pos in c1.execute(cblu):
        calls[ii][pos] = 1
    trace(1, 'cblu found at {:.2f} seconds'.format(time.time() - t0))
    for ii,pos in c1.execute(cbu):
        calls[ii][pos] = 2
    trace(1, 'cbu found at {:.2f} seconds'.format(time.time() - t0))
    for ii,pos in c1.execute(cbl):
        calls[ii][pos] = 3
    trace(1, 'cbl found at {:.2f} seconds'.format(time.time() - t0))

    # get list of ids for the kits, in their original listing order
    iids = []
    files = c1.execute('select id,name,seq from files where kit=1 order by 3')
    for ii,fname,seq in files:
        iids.append(ii)

    # track the pos,ref,alt,id of the variants
    vrai_vector = []
    trace(3,'vcounts: {}'.format(vcounts))
    for vid in vcounts:
        c1 = c1.execute('select pos,ref,alt from variants where id=?', (vid,))
        pos,ref,alt = next(c1)
        vrai_vector.append((pos,ref,alt,vid))
        trace(3,'{}'.format(vrai_vector[-1]))
    vrai_vector.sort()
    vra_vector = [(a,b,c) for (a,b,c,d) in vrai_vector]
    v_vector = [a for (a,b,c,d) in vrai_vector]
    vid_vector = [d for (a,b,c,d) in vrai_vector]

    # all of the calls across the kits
    kit_calls = defaultdict(dict)
    kcalls = c1.execute('select id,vid from vcfcalls')
    for ID,vid in kcalls:
        kit_calls[ID][vid] = None
    trace(1, 'raw variants done at {:.2f} seconds'.format(time.time() - t0))

    # find if variant was covered (note cbu/cbl already handled above)
    ist = 0
    for ii in iids:
        ist += 1
        ranges = c1.execute('select minaddr,maxaddr from bed where id=? order by minaddr', (ii,))
        kit_ranges = []
        for minaddr,maxaddr in ranges:
            kit_ranges.append((minaddr,maxaddr))

        coverage_vector = in_range(v_vector, kit_ranges)
        for iv,v in enumerate(v_vector): # index,address
            if vid_vector[iv] in kit_calls[ii]:
                continue # will write in output below
            if not coverage_vector[iv] and v not in calls[ii]:
                calls[ii][v] = 4 # ;nc
        if (ist % 250) == 249:
            trace(2, '250 people done at {:.2f} seconds'.format(time.time() - t0))

    trace(1, '{} rows done at {:.2f} seconds'.format(len(vra_vector), time.time() - t0))

    # rejects across all of the kits
    rejects = defaultdict(dict)
    rej = dbcurs.execute('select id,vid from vcfrej')
    for iid,vid in rej:
        rejects[iid][vid] = None

    # write the output
    out = ''
    for iiv, vid in enumerate(vid_vector):
        vname = '{}'.format(v_vector[iiv])
        if vid in ref_swaps and vcounts[vid] > 0:
            nvar = 0
            for ii in iids:
                # a kit that was not called for this variant...reverse
                if vid not in kit_calls[ii] and ((v_vector[iiv] not in calls[ii]) or (calls[ii][v_vector[iiv]] != 4)):
                    if vid not in rejects[ii]:
                        nvar += 1
            ref,alt = ref_swaps[vid]
            trace(2,'swapping variant id={}; nvar={}'.format(vid,nvar))
        else:
            nvar = vcounts[vid]
            xx,ref,alt = vra_vector[iiv]
        if nvar == 0 and vid not in inserted:
            continue # skip over swap variants, but keep inserted variants
        ref = ref.replace(',', '|') # cleanse any commas
        alt = alt.replace(',', '|')
        out += vname+',,'+ref+','+alt+','
        if vid in snps:
            out += 'SNP,'
        elif vid in indels:
            out += 'Indel,'
        out += '{},,,,,,,,,,,'.format(nvar)
        for ii in iids:
            out += ','
            if vid in ref_swaps and vcounts[vid] > 0:
                if vid in kit_calls[ii]:
                    tcall = ''
                elif (v_vector[iiv] not in calls[ii] or calls[ii][v_vector[iiv]] != 4):
                    if vid not in rejects[ii]:
                        tcall = vname+'.'+ref+'.'+alt
                    else:
                        tcall = ';rej' # mark rejects differently
                        tcall = '' # compatibility with existing script
                else:
                    tcall = ''
            elif vid in kit_calls[ii]:
                tcall = vname+'.'+ref+'.'+alt
            else:
                tcall = ''
            out += tcall
            if v_vector[iiv] in calls[ii]:
                out += call_strings[calls[ii][v_vector[iiv]]]
        print(out)
        out = ''

    trace (1, 'all done at {:.2f} seconds'.format(time.time() - t0))



# STATS1 in redux.bash
if stats1:
    c1 = dbconn.cursor()
    c = c1.execute('select id, coverage1, coverage2, nranges from bedstats')
    for row in c:
        print(row[1], row[2], row[3])



# STATS2 in redux.bash
if stats2:
    c1 = dbconn.cursor()
    c = c1.execute('select id, ny, nv, ns, ni from vcfstats')
    for row in c:
        print(row[1], row[2], row[3], 0,0,0, row[4], 0,0)



# list out which files are stored, in sort order
if listfiles:
    c1 = dbconn.cursor()
    c = c1.execute('select name,seq from files where kit=1 order by 2')
    for row in c:
        print(row[0])



# list out full path bed files, in sort order
if listbed:
    c1 = dbconn.cursor()
    c = c1.execute('select name,seq from files where kit=1 order by 2')
    for row in c:
        print(os.path.join(UNZIPDIR,row[0])+'.bed')



# update the snp defs and names from hg19
if updatesnps:
    with open('snps_hg19.csv') as snpfile:
        c1 = dbconn.cursor()
        reader = csv.DictReader(snpfile)
        snps = defaultdict(list)
        for row in reader:
            snpdef = (row['start'], row['allele_anc'], row['allele_der'])
            snps[snpdef].append(row['Name'])
        for snp in snps:
            c1.execute('insert or ignore into variants(pos,ref,alt) values(?,?,?)', snp)
            c1 = c1.execute('select id from variants where pos=? and ref=? and alt=?', snp)
            vid = c1.fetchone()[0]
            for n in snps[snp]:
                c1.execute('insert into snpnames values(?,?)', (vid, n))



# print info about a snp by name or addr
if querysnp:
    trace(2, 'looking for details about {}...'.format(querysnp))
    c1 = dbconn.cursor()
    ids = set()
    c1 = c1.execute('select * from variants where pos=?', (querysnp,))
    for row in c1:
        ids.add(row[0])
    # wildcard version - usually not useful
    #c1 = c1.execute('select * from snpnames where name like ?', (querysnp+'%',))
    c1 = c1.execute('select id from snpnames where name=?', (querysnp.upper(),))
    for row in c1:
        ids.add(row[0])

    poss = set()
    nams = set()
    for id in ids:
        c1 = c1.execute('select distinct pos from variants where id=?', (id,))
        for row in c1:
            poss.add('{}'.format(row[0]))
        c1 = c1.execute('select distinct name from snpnames where id=?', (id,))
        for row in c1:
            nams.add(row[0])
    for pos in poss:
        c1 = c1.execute('select pos,ref,alt,id from variants where pos=?', (pos,))
        for row in c1:
            print('{:>8}.{}.{} - {} (id={})'.format(row[0], row[1], row[2], '/'.join(nams),row[3]))

    if not ids:
        print('No data on', querysnp)

    kitids = set()
    if vars(namespace)['kits']:
        kitids = set()
        for iid in ids:
            c1 = c1.execute('select distinct id from vcfcalls where vid=?', (iid,))
            for row in c1:
                kitids.add(row[0])
        print('Found in', len(kitids), 'kits')
        for iid in kitids:
            c1 = c1.execute('select name from files where id=?', (iid,))
            print(c1.fetchone()[0])

    if vars(namespace)['implications']:
        if poss:
            print('\n--- implications.txt ---')
            io = subprocess.run(['egrep', '-iw', '|'.join(poss), 'implications.txt'],
                            stdout=subprocess.PIPE)
            print(io.stdout.decode('utf-8'))

    if vars(namespace)['tree']:
        if poss:
            print('--- tree.csv ---')
            io = subprocess.run(['egrep', '-iw', '|'.join(poss), 'tree.csv'],
                                stdout=subprocess.PIPE)
            print(io.stdout.decode('utf-8'))

    if vars(namespace)['badlist']:
        if poss:
            print('--- badlist.txt ---')
            io = subprocess.run(['egrep', '-iw', '|'.join(poss), 'badlist.txt'],
                                stdout=subprocess.PIPE)
            print(io.stdout.decode('utf-8'))

            print('--- recurrencies.txt ---')
            io = subprocess.run(['egrep', '-iw', '|'.join(poss), 'recurrencies.txt'],
                                stdout=subprocess.PIPE)
            print(io.stdout.decode('utf-8'))

#incomplete
if mergeup:
    c1 = dbconn.cursor()
    c2 = dbconn.cursor()
    c3 = dbconn.cursor()
    # populate tree from the rows in tree csv
    # store the variants associated with the row (column 5, semicolon sep'd)
    c1.execute('CREATE TEMPORARY TABLE tree(e INTEGER PRIMARY KEY, a INTEGER, b INTEGER, c INTEGER, d INTEGER, f INTEGER, g TEXT, h TEXT, i INTEGER)')
    c1.execute('CREATE TEMPORARY TABLE treevars(ID INTEGER REFERENCES tree(e), snpid INTEGER REFERENCES variants(id), unique(id,snpid))')
    with open('tree.csv') as csvfile:
        fieldnames=['a','b','c','d','e','f','g','h','i']
        reader = csv.DictReader(csvfile, fieldnames=fieldnames)
        for row in reader:
            tup=(row['a'],row['b'],row['c'],row['d'],
                     row['f'],row['g'],row['h'],row['i'])
            c1.execute('insert into tree(a,b,c,d,f,g,h,i) values (?,?,?,?,?,?,?,?)', tup)
            treeid=c1.lastrowid
            snplist = row['e'].split(';')

            for snp in snplist:
                if '/' in snp:
                    poss = set()
                    for sname in snp.split('/'):
                        c1 = c1.execute('select id from snpnames where name=?', (sname,))
                        r1 = next(c1)
                        if r1:
                            c1.execute('insert or ignore into  treevars values(?,?)', (treeid,snpid))
                else:
                    s = c1.execute('select id from variants where pos=?', (snp,))
                    try:
                        snpid = s.fetchone()[0]
                    except:
                        c1.execute('insert into variants(pos) values(?)', (snp,))
                        snpid = c1.lastrowid
                    c1.execute('insert or ignore into treevars values(?,?)', (treeid,snpid))

    # Each row in the tree carries a range [minkit, maxkit] (col3 & col4)
    # Those are the column indices in the variant-not-shared.txt array.
    # Columns represent people, so it's a range of people in the clade.
    # For a sub-node to be merged up, it is necessary that all of the
    # people grouped in the sub-node's range are also in the parent.
    # So we have to check the parent range of people against all of the snps
    # that define the child clade, and if any such snp is NOT
    # called for some portion of the people in the parent range,
    # it precludes merging.
    # example:
    # kits k2,k3,k4,k5 are part of a parent clade
    # kits k4,k5 have snps S1,S2 that k2,k3 don't have
    # snps S1,S2 define a child clade that can't be merged up to the parent

    # in this table, let's make person id = column number
    c1.execute('CREATE TEMPORARY TABLE vars(snp INTEGER REFERENCES variants(id), pid INTEGER)')
    tups = []
    with open('variant-not-shared.txt') as varfile:
        c = csv.reader(varfile)
        for row in c:
            rowsnp = row[0]
            try:
                rowsnpid = c2.execute('select id from variants where pos=?', (rowsnp,)).fetchone()[0]
            except:
                # insert variant since we don't have it yet
                trace(2, 'pos = {} added'.format(row[0]))
                c2.execute('insert into variants(pos) values(?)', (rowsnp,))
                rowsnpid = c2.lastrowid
            for col in range(17,len(row)):
                if row[col]:
                    tups.append((rowsnpid, col+1))
    c1.executemany('insert into vars values(?,?)', tups) 
    # tables populated

    c1 = c1.execute('select * from tree where g!="0"')
    for row in c1:
        par = os.path.splitext(row[6])[0]
        prow = next(c2.execute('select e,c,d from tree where g=?', (par,)))
        minidx = prow[1]
        maxidx = prow[2]
        trace(3,'id:{}; max/min:{}/{}'.format(prow[0], minidx, maxidx))
        c2 = c2.execute('select snpid from treevars where id=?', (row[0],))
        for r1 in c2:
            for kit in range(minidx,maxidx+1):
                trace(3,'checking kit {} for snpid {}'.format(kit,r1[0]))
                hasit = c3.execute('select id from vcfcalls where vid=? and id=(select id from files where seq=?)', (r1[0], kit))
                try:
                    snpid = next(hasit)
                    continue
                except:
                    break # cannot be merged up
            else:
                continue
            break
        else:
            print(row[6], 'into its parent')
            continue # on to the next row
            

dbconn.commit()
dbcurs.close()
trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))

