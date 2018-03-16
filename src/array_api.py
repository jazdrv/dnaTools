#!/usr/bin/env python3
# coding: utf-8

# Contributors: Jef Treece, Harald Alvestrand, Zak Jones, Iain McDonald
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

from collections import defaultdict
import sys, yaml, time, os
from lib import Trace, unpack_call

# read the config file
sys.path.insert(0, os.environ['REDUX_PATH'])
REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
config = yaml.load(open(REDUX_CONF))

trace = Trace(config['verbosity'])

# Procedure: get_kit_coverage
# Purpose:
#   return a) sum of BED ranges for kit (how many bases are covered by the
#   test) and b) sum of coverage that is in the age BED range. These stats are
#   needed for age calculations.
# Input:
#   dbo, a database object
#   pid, a person ID to be updated
# Returns:
#   total sum of coverage (sum of bed ranges)
#   total sum of coverage that is within the agebed ranges
# Info:
#   This is not done in the "brute force" way because it can be compute
#   intensive to search the list for every range.
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

# Procedure: in_range
# Purpose: determine if a set of positions is contained in a set of ranges
# Input:
#   v_vect, a sorted vector of positions
#   ranges, a range vector sorted on minaddr (minaddr,maxaddr)
# Returns:
#   a vector of True values if v contained in a range, else False
# Info:
#   does not consider the endpoints; returns True if minaddr < pos < maxaddr
def in_range(v_vect, ranges, spans):
    trace(4,'in_range')
    c_vect = []
    ii = 0
    nmax = len(ranges)
    for iv,v in enumerate(v_vect):
        while ii < nmax and ranges[ii][0] < v:
            ii += 1
        if ii > 0 and v > ranges[ii-1][0] and v < ranges[ii-1][1]:
            # note <= ranges[1] means we don't treat cbh any differently
            if spans and v+spans[iv] <= ranges[ii-1][1]:
                # span fits within the range
                c_vect.append(True)
                trace(5, '{}-{}:cov ({},{})'.format(iv,v,ranges[ii-1][0],ranges[ii-1][1]))
            elif spans:
                # span exceeded the upper end of the range
                c_vect.append(False)
                trace(5, '{}-{}:nc ({},{})'.format(iv,v,ranges[ii-1][0],ranges[ii-1][1]))
            else:
                # no span - SNP fits within range
                c_vect.append(True)
                trace(5, '{}-{}:snp cov ({},{})'.format(iv,v,ranges[ii-1][0],ranges[ii-1][1]))
        else:
            c_vect.append(False)
            trace(5, '{}-{}:not cov ({},{})'.format(iv,v,ranges[ii-1][0],ranges[ii-1][1]))
    if len(c_vect) != len(v_vect):
        raise ValueError
    return c_vect


# Procedure: get_variant_array
# Purpose: build a dictionary of people-variants
# Input:
#   db, a database object
#   ppl, a 1-d array of dataset IDs we're interested in
# Returns:
#   the variant array, a 2d dictionary, True if person has the variant
#   list of people represented
#   list of variants represented
# Info:
#   All entries in vcfcalls are returned; not limited to repeated calls.  Call
#   get_kit_coverages in conjunction with this to get the coverage associated
#   with the list of variants
def get_variant_array(db, ppl):
    # FIXME: handle refpos
    # FIXME: handle multiple calls at a given pos for a given person?
    # FIXME: downstream analysis may need more than True/False
    # FIXME: merge identical indels here?
    c1 = db.dc
    c1.execute('drop table if exists tmpt')
    c1.execute('create temporary table tmpt(id integer)')
    c1.executemany('insert into tmpt values(?)', [(p,) for p in ppl])
    c1 = c1.execute('''select c.vid,c.pid,c.callinfo from vcfcalls c
                     inner join tmpt t on c.pid=t.id''')
    # build a 2d array
    arr = {}
    var = set()
    for pp in ppl:
        arr[pp] = defaultdict()

    for (v,p,c) in c1:
        passfail, q1, q2, numcalls, passrate = unpack_call(c)
        if passfail and numcalls > 1:
            arr[p][v] = True
        var.add(v)

    c1.execute('drop table tmpt')

    return arr, ppl, list(var)

# Procedure: get_variant_csv
# Purpose: create a csv file of 2d person x variant array
# Inputs:
#   db, a database object
#   ppl, a 1-d array of dataset IDs we're interested in
# Returns:
#   a csv file as a string
# Info:
#   Materializing this array is impractical if ppl is very large.  This
#   procedure is useful for human-readable debugging.  Cells in the csv have +
#   if the variant is called and ';nc' if it's not covered by a BED range.
def get_variant_csv(db, ppl):
    trace(2, 'get_variant_csv at {}'.format(time.clock()))
    # people of interest
    dc = db.cursor()
    dc.execute('drop table if exists tmpp')
    dc.execute('create temporary table tmpp(pid integer)')
    dc.executemany('insert into tmpp(pid) values(?)', [(p,) for p in ppl])

    # get variants in ascending order of number of calls across kits
    vc = db.cursor()
    vc.execute('''select v.id,v.anc,v.der,count(c.pid)
                  from variants v
                  inner join vcfcalls c on c.vid = v.id
                  inner join tmpp t on t.pid = c.pid
                  group by 1,2,3 order by 4''')

    # variant list is used more than once, so make a list
    # keep only variants that have at least two call entries across all ppl
    # this may need to change when considering singletons
    vl = list([s[0] for s in vc if s[3] > 1])

    cdict = get_kit_coverages(db, ppl, vl)
    tc = db.cursor()
    tc.execute('drop table if exists tmpv')
    tc.execute('create temporary table tmpv(vid integer)')
    tc.executemany('insert into tmpv(vid) values(?)', [(v,) for v in vl])
    tc.execute('''select v.id,v.pos,a.allele,b.allele from variants v
                  inner join alleles a on a.id=v.anc
                  inner join alleles b on b.id=v.der
                  inner join tmpv t on t.vid=v.id''')
    vlist=list([(v[0],'{}.{}.{}'.format(v[1],v[2],v[3])) for v in tc])

    # get all calls appearing in vcfcalls
    arr, pl, vlfull = get_variant_array(db, ppl)
    tc.execute('''select d.id,d.kitid from dataset d
                  inner join tmpp t on t.pid=d.id''')
    klist = list([(k[0],k[1]) for k in tc])

    trace(2, 'write csv file at {}'.format(time.clock()))
    out = ',' + ','.join(['{}'.format(v[1]) for v in klist]) + '\n'
    for vid,vtext in vlist:
        out += vtext
        out += ','
        for pid,kid in klist:
            # if there's an entry, it means not covered
            try:
                coverage = cdict[pid][vid]
                ncstring = ';nc'
            except:
                # no entry in the array means covered
                ncstring = ''
            try:
                # entry in the array means called
                callinfo = arr[pid][vid]
                callstring = '+'
            except:
                callstring = ''

            out += '{}{},'.format(callstring, ncstring)
        out += '\n'
    return out

# Procedure: get_dna_ids
# Purpose: get the list of populated (ones that have vcfcalls entries) DNAIDs
# Input: db, a database object
def get_dna_ids(db):
    return [x[0] for x in db.dc.execute('select distinct(pID) from vcfcalls')]

# Procedure: get_call_coverage
# Purpose: check coverage for a vector of variants for a given DNA kit
# Input:
#   dbo, a database object
#   pid, a person ID
#   vids, a vector of variant ids to be checked
#   spans, optional vector of spans associated with vids; a SNP span is 1
# Info:
#   Check kit coverage for a vector of variants; return a coverage vector that
#   indicates if that person has a BED range such that the call is in a
#   range. If spans is passed, it corresponds to the maximum length affected by
#   the variant in vids: 1 for SNPs, max length of (ref,alt) for indels
def get_call_coverage(dbo, pid, vids, spans=None):
    # FIXME currently only returns True/False, doesn't handle range ends
    # FIXME maybe more efficient to pass positions instead of ids
    # FIXME maybe calculate spans from variants instead of passing it in
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
    xl = list([(v[0],v[1]) for v in calls])
    iv = [v[0] for v in xl]
    pv = [v[1] for v in xl]
    trace(500, '{} calls: {}...'.format(len(pv), pv[:20]))
    # form a list of ranges of interest
    rc = dbo.cursor()
    ranges = rc.execute('''select minaddr,maxaddr from bedranges r
                           inner join bed b on b.bID=r.id
                           where b.pid=?
                           order by 1''', (pid,))
    rv = [v for v in ranges]
    if len(rv) == 0:
        return []
    trace(500, '{} ranges: {}...'.format(len(rv), rv[:20]))
    # calculate coverage vector
    coverage = in_range(pv, rv, spans)
    dc.close()
    rc.close()
    return iv, coverage

# Procedure: get_kit_coverages
# Purpose: find coverage for a list of kits and a list of variants
# Input:
#   db, a database object
#   pids, a vector of person ids
#   vids, a vector of variant ids
# Returns:
#   2d dictionary indexed by pid,vid with False if vid not covered
# Info:
#   determines indels and snps in the list of variants and calls
#   get_call_coverage for each person in the list
def get_kit_coverages(db, pids, vids):

    # these are the variants to inspect
    c = db.cursor()
    c.execute('drop table if exists tmpt')
    c.execute('create table tmpt(vid integer)')
    c.executemany('insert into tmpt(vid) values(?)', [(v,) for v in vids])

    c = db.cursor()
    c.execute('''select distinct v.id,max(length(a.allele),length(b.allele))
                 from variants v
                 inner join alleles a on a.id=v.anc
                 inner join alleles b on b.id=v.der
                 inner join tmpt c on c.vid=v.id
                 where length(a.allele)>1 or length(b.allele)>1''')
    indels = [t for t in c]
    indel_ids = list([t[0] for t in indels])
    spans = list([int(t[1]) for t in indels])

    c = db.cursor()
    c.execute('''select distinct v.id
                 from variants v
                 inner join alleles a on a.id=v.anc
                 inner join alleles b on b.id=v.der
                 inner join tmpt c on c.vid=v.id
                 where length(a.allele)=1 and length(b.allele)=1''')
    snps = [t for t in c]
    snp_ids = list([t[0] for t in snps])

    cdict = {}
    trace(1, 'calculate coverages')
    for pid in pids:
        cdict[pid] = defaultdict()
    for pid in pids:
        # get indel coverage for a kit
        trace(3, 'indels for kit {} at {}...'.format(pid,time.clock()))
        trace(4, 'get_call_coverage(db, {}, {}, {})'.format(pid,[(i[0],i[1]) for i in enumerate(indel_ids)][:50],[(i[0],i[1]) for i in enumerate(spans)][:50]))
        iv,cv = get_call_coverage(db, pid, indel_ids, spans)
        trace(5, 'indels:{}..., coverage:{}...'.format([(i[0],i[1]) for i in enumerate(iv)][:50], [(i[0],i[1]) for i in enumerate(cv)][:50]))
        # store "not-covered" since it's sparse
        for cov,vid in zip(cv, iv):
            if not cov:
                cdict[pid][vid] = False
        # get snp coverage for a kit
        trace(3, 'snps for kit {} at {}...'.format(pid,time.clock()))
        trace(4, 'get_call_coverage(db, {}, {})'.format(pid,snp_ids))
        iv,cv = get_call_coverage(db, pid, snp_ids)
        trace(5, 'snps:{}..., coverage:{}...'.format(iv[:5], cv[:5]))
        # store "not-covered" since it's sparse
        for cov,vid in zip(cv, iv):
            if not cov:
                cdict[pid][vid] = False
    return cdict

# test framework
if __name__=='__main__':
    trace(0, 'test message should display to stdout')
