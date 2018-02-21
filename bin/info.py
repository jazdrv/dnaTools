#!/usr/bin/env python3
"""
  Purpose:
    report information from the redux database

  Usage:
    -s <snp>  show information about a snp by name or position

  Copyright:
    For free distribution under the terms of the
    GNU General Public License, version 3 (29 June 2007)
    https://www.gnu.org/licenses/gpl.html

  Jef Treece, 21 Feb 2018
"""

import sqlite3, yaml, sys, time, argparse

t0 = time.time()

config = yaml.load(open('config.yaml'))

# diagnostics
def trace (level, msg, stream=sys.stderr):
    if level <= config['verbosity']:
        if level == 0:
            print(msg)
        else:
            print(msg, file=stream)
            stream.flush()

# affects diagnostic messages, which always go to stderr
# this is only the default; increase with -v flag(s)
DEBUG=config['verbosity']

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--snp', nargs=1)
args = parser.parse_args()


# create a sqlite database, cursor, and tables
dbconn = sqlite3.connect(config['DB_FILE'])
dbcurs = dbconn.cursor()


# STATS1 in redux.bash
def stats1():
    c1 = dbconn.cursor()
    c = c1.execute('select id, coverage1, coverage2, nranges from bedstats')
    for row in c:
        print(row[1], row[2], row[3])


# STATS2 in redux.bash
def stats2():
    c1 = dbconn.cursor()
    c = c1.execute('select id, ny, nv, ns, ni from vcfstats')
    for row in c:
        print(row[1], row[2], row[3], 0,0,0, row[4], 0,0)


# list out which files are stored, in sort order
def listfiles():
    c1 = dbconn.cursor()
    c = c1.execute('select name,seq from files where kit=1 order by 2')
    for row in c:
        print(row[0])


# list out full path bed files, in sort order
def listbed():
    c1 = dbconn.cursor()
    c = c1.execute('select name,seq from files where kit=1 order by 2')
    for row in c:
        print(os.path.join(UNZIPDIR,row[0])+'.bed')


# print info about a snp by name or addr
def querysnp(snp):
    trace(2, 'looking for details about {}...'.format(snp))
    c1 = dbconn.cursor()
    ids = set()
    c1 = c1.execute('select id from variants where pos=?', (snp,))
    for row in c1:
        ids.add(row[0])
    c1 = c1.execute('select vid from snpnames where snpname=?', (snp.upper(),))
    for row in c1:
        ids.add(row[0])

    poss = set()
    nams = set()
    for id in ids:
        c1 = c1.execute('select distinct pos,buildid from variants where id=?',
                            (id,))
        for row in c1:
            poss.add('{}'.format(row[0]))
        c1 = c1.execute('select distinct snpname from snpnames where vid=?',
                            (id,))
        for row in c1:
            nams.add(row[0])
    for pos in poss:
        c1 = c1.execute('''select bld.buildnm, v.pos, a.allele as anc,
                           b.allele as der, v.id
                           from variants v
                           inner join alleles a on a.id=v.anc
                           inner join alleles b on b.id=v.der
                           inner join build bld on bld.id=v.buildid
                           where pos=?''', (pos,))
        for row in c1:
            print('{} {:>8}.{}.{} - {} (id={})'.format(row[0], row[1], row[2],
                                                       row[3],
                                                       '/'.join(nams),
                                                       row[4]))
    if not ids:
        print('No data on {}'.format(snp))

    kitids = set()
    for iid in ids:
        c1 = c1.execute('select distinct pid from vcfcalls where vid=?', (iid,))
        for row in c1:
            kitids.add(row[0])
    print('Found in {} kits'.format(len(kitids)))
    for iid in kitids:
        c1 = c1.execute('''select d.kitid,b.buildnm from dataset d
                           inner join build b on b.id=d.buildid
                           where d.id=?''', (iid,))
        kitid, buildnm = c1.fetchone()
        print('{} {}'.format(kitid,buildnm))


if args.snp:
    querysnp(args.snp[0])

dbconn.commit()
dbcurs.close()
trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))

