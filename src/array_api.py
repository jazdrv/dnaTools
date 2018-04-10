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

# values that may be returned by in_range
RANGE_VALS = 'nc', 'cbl', 'cbh', 'cblh', 'cov'
RANGE_NC = 0
RANGE_CBL = 1
RANGE_CBH = 2
RANGE_CBLH = 3
RANGE_COV = 4

# Procedure: in_range
# Purpose: determine if a set of positions is contained in a set of ranges
# Input:
#   v_vect, a sorted vector of positions
#   ranges, a range vector sorted on minaddr (minaddr,maxaddr)
#   spans, a vector of lengths associates with v_vect or None
# Returns:
#   a vector of values:
#     0: if v is not contained in a range
#     1: if v is equal to the lower bound of a range
#     2: if v is equal to the upper bound of a range
#     3: if v equals both upper and lower bound of a range (of one)
# Info:
#   Spans are for indels, which span a range of positions.
#   FTDNA ranges are zero-based.
#   For example, the range (1,2) refers to position 2 only, and (5,10) refers
#   to position 6-10 inclusive. In this case, a SNP at position 6 would be cbl.
def in_range(v_vect, ranges, spans):
    trace(4,'in_range')
    c_vect = []
    ii = 0
    nmax = len(ranges)
    for iv,v in enumerate(v_vect):
        while ii < nmax and (ranges[ii][0]+1) <= v:
            ii += 1
        if not spans: # snps
            if v == ranges[ii-1][1] and v == (ranges[ii-1][0] + 1):
                c_vect.append(RANGE_CBLH)
                trace(25, '{}:cblh ({},{})'.format(v,ranges[ii-1][0],
                                                       ranges[ii-1][1]))
            elif v == ranges[ii-1][1]:
                c_vect.append(RANGE_CBH)
                trace(25, '{}:cbh ({},{})'.format(v,ranges[ii-1][0],
                                                      ranges[ii-1][1]))
            elif v == (ranges[ii-1][0]+1):
                c_vect.append(RANGE_CBL)
                trace(25, '{}:cbl ({},{})'.format(v,ranges[ii-1][0],
                                                      ranges[ii-1][1]))
            elif ii > 0 and v > ranges[ii-1][0] and v < ranges[ii-1][1]:
                c_vect.append(RANGE_COV)
                trace(25, '{}:snp cov ({},{})'.format(v,ranges[ii-1][0],
                                                          ranges[ii-1][1]))
            else:
                c_vect.append(RANGE_NC)
                trace(25, '{}:snp nc ({},{})'.format(v,ranges[ii-1][0],
                                                         ranges[ii-1][1]))
        else: # spans
            if v == (ranges[ii-1][0]+1) and v+spans[iv]-1 == ranges[ii-1][1]:
                c_vect.append(RANGE_CBLH)
                trace(25, '{}-{}:cblhspan ({},{})'.format(v,v+spans[iv]-1,
                                            ranges[ii-1][0], ranges[ii-1][1]))
            elif v == (ranges[ii-1][0]+1) and v+spans[iv]-1 < ranges[ii-1][1]:
                c_vect.append(RANGE_CBL)
                trace(25, '{}-{}:cblspan ({},{})'.format(v,v+spans[iv]-1,
                                            ranges[ii-1][0], ranges[ii-1][1]))
            elif v+spans[iv]-1 == ranges[ii-1][1]:
                c_vect.append(RANGE_CBH)
                trace(25, '{}-{}:cbhspan ({},{})'.format(v,v+spans[iv]-1,
                                            ranges[ii-1][0], ranges[ii-1][1]))
            elif v+spans[iv]-1 < ranges[ii-1][1]:
                # span fits within the range
                c_vect.append(RANGE_COV)
                trace(25, '{}-{}:covspan ({},{})'.format(v,v+spans[iv]-1,
                                            ranges[ii-1][0], ranges[ii-1][1]))
            else:
                # span exceeded the upper end of the range
                c_vect.append(RANGE_NC)
                trace(25, '{}-{}:ncspan ({},{})'.format(v,v+spans[iv]-1,
                                            ranges[ii-1][0], ranges[ii-1][1]))
    if len(c_vect) != len(v_vect):
        raise ValueError
    return c_vect


# Procedure: get_variant_defs
# Purpose: get the ref and alt ids for a list of variants
# Input:
#   a db instance
#   a dict of dnaid:kitid
# Returns:
#   (vid, pos, anc, der) for vid in variants
def get_variant_defs(db, vids):
    dc = db.cursor()
    rval = {}
    for v in vids:
        dc.execute('''select v.id,v.pos,aa.allele,ab.allele from variants v
                    inner join alleles aa on v.anc=aa.id
                    inner join alleles ab on v.der=ab.id
                    where v.id=?''', (v,))
        tup = dc.fetchone()
        rval[tup[0]] = (tup[1], tup[2], tup[3])
    return rval


# Procedure: get_variant_defs
# Purpose: get the ref and alt ids for a list of variants
# Input:
#   a db instance
#   a dict of dnaid:kitid
# Returns:
#   (vid, pos, anc, der) for vid in variants
def get_variant_snpnames(db, vids):
    dc = db.cursor()
    rval = {}
    for v in vids:
        dc.execute('''select v.id,s.snpname from variants v
                    inner join snpnames s on s.vid=v.id
                    where v.id=?''', (v,))
        sn = '/'.join([s[1] for s in dc])
        # note sn might be empty, one, or multiple names, e.g. 'M2986/Z4303'
        rval[v] = sn
    return rval


# Procedure: get_kit_ids
# Purpose: get the lab-assigned kitid corresponding to dnaid
# Input:
#   a db instance
#   a dict of dnaid:kitid
# Returns:
#   (dnaid, kitid) for dnaid in pids
def get_kit_ids(db, pids):
    dc = db.cursor()
    rval = {}
    for p in pids:
        dc.execute('select dnaid,kitid from dataset where dnaid=?', (p,))
        pid,kit = dc.fetchone()
        rval[pid] = kit
    return rval


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
#   with the list of variants returned.
#   
#   Note that this can be called without the ppl argument, in which case ppl is
#   formed by querying the analysis_kits table
def get_variant_array(db, ppl=None):
    # FIXME: handle multiple calls at a given pos for a given person?
    # FIXME: downstream analysis may need additional info in return tuple
    # FIXME: merge identical indels here?
    c1 = db.dc
    c1.execute('drop table if exists tmpt')
    c1.execute('create temporary table tmpt(id integer)')
    if ppl:
        c1.executemany('insert into tmpt values(?)', [(p,) for p in ppl])
    else:
        c1.execute('insert into tmpt select * from analysis_kits')
    c1 = c1.execute('''select c.vid,c.pid,c.callinfo from vcfcalls c
                     inner join tmpt t on c.pid=t.id''')

    # refpos processing
    c2 = db.db.cursor()
    c2.execute('''select distinct r.vid, rv.id from refpos r
                  inner join variants v on v.id=r.vid
                  inner join variants rv on rv.anc=v.der and
                     rv.der=v.anc and rv.pos=v.pos and rv.buildid=v.buildid
               ''')
    refpos = dict([(rp[0],rp[1]) for rp in c2])
    c2.close()

    # build a 2d array
    arr = {}
    var = set()
    for pp in ppl:
        arr[pp] = defaultdict()

    for (v,p,c) in c1:
        passfail, gt, q1, q2, numcalls, passrate = unpack_call(c)
        # handle refpos variants
        try:
            # swap variant and call's genotype
            v = refpos[v]
            gt = {0:1, 1:0, 2:2, 3:3}[gt]
        except:
            pass
        if passfail and numcalls > 1:
            arr[p][v] = (True, gt)
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
    tc.execute('''select d.dnaid,d.kitid from dataset d
                  inner join tmpp t on t.pid=d.dnaid''')
    klist = list([(k[0],k[1]) for k in tc])

    trace(2, 'write csv file at {}'.format(time.clock()))
    out = ',' + ','.join(['{}'.format(v[1]) for v in klist]) + '\n'
    for vid,vtext in vlist:
        out += vtext
        out += ','
        for pid,kid in klist:
            try:
                # indicate what coverage (blank means covered)
                if cdict[pid][vid] < RANGE_COV:
                    ncstring = RANGE_VALS[cdict[pid][vid]]
                else:
                    ncstring = ''
            except:
                ncstring = ''
            try:
                # entry in the array means called
                callinfo = arr[pid][vid]
                callstring = '+'
            except:
                callstring = ''
            if ncstring and callstring:
                out += ';'.join([callstring, ncstring])
            else:
                out += (callstring + ncstring)
            out += ','
        out = out[:-1] + '\n'
    return out

# Procedure: get_dna_ids
# Purpose: get the list of populated (ones that have vcfcalls entries) DNAIDs
# Returns: a vector of pIDs that have calls
# Input: db, a database object
def get_dna_ids(db):
    return [x[0] for x in db.dc.execute('select distinct(pID) from vcfcalls')]

# Procedure: get_analysis_ids
# Purpose: get the list of populated DNAIDs for analysis
# Input: db, a database object
# Returns: a vector of pIDs in analysis_kits that also have data (intersection)
# Info:
#   We may have more kits loaded than we want to analyze. The table
#   analysis_kits has a list of IDs we want to analyze. This permits loading
#   arbitrary kits into the database, while keeping the analysis focused on a
#   subset of kits.
def get_analysis_ids(db):
    all_ids = set(get_dna_ids(db))
    trace(5, 'all ids: {}'.format(all_ids))
    ids = set([x[0] for x in db.dc.execute('select pID from analysis_kits')])
    trace(5, 'analysis ids: {}'.format(ids))
    return list(ids.intersection(all_ids))

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
    # get list of variants of interest ordered by position
    calls = dc.execute('''select v.id, v.pos from variants v
                          inner join tmpt t on t.vid=v.id
                          order by 2''')
    # form a list of positions of interest
    xl = list([(v[0],v[1]) for v in calls])
    iv = [v[0] for v in xl]
    pv = [v[1] for v in xl]
    trace(500, '{} calls: {}...'.format(len(pv), pv[:20]))
    # get the ordered list of ranges of interest
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
        # store coverage information
        for cov,vid in zip(cv, iv):
            if cov != RANGE_COV:
                cdict[pid][vid] = cov
        # get snp coverage for a kit
        trace(3, 'snps for kit {} at {}...'.format(pid,time.clock()))
        trace(4, 'get_call_coverage(db, {}, {})'.format(pid,snp_ids))
        iv,cv = get_call_coverage(db, pid, snp_ids)
        trace(5, 'snps:{}..., coverage:{}...'.format(iv[:5], cv[:5]))
        # store "not-covered" since it's sparse
        for cov,vid in zip(cv, iv):
            if cov != RANGE_COV:
                cdict[pid][vid] = cov
    return cdict

# test framework
if __name__=='__main__':
    from db import DB
    from subprocess import call
    t0 = time.time()
    # smoke test DB __init__
    db = DB(drop=False)
    # smoke test trace() and get_analysis_ids()
    trace(0, 'test message should display to stdout')
    ids = get_analysis_ids(db)
    trace(0, 'analysis_ids: {} (empty if kits are not loaded)'.format(ids))
    for id in ids:
        call(['../bin/info.py', '-k', '{}'.format(id)])
    # smoke test in_range()
    v_vect = [1, 1, 5, 30, 35, 40, 42, 47, 52]
    ranges = [(0, 20), (30, 40), (41, 42), (45, 50)]
    # expected: [1, 1, 4, 0, 4, 2, 3, 4, 0]
    trace(0, '{}'.format(in_range(v_vect,ranges,None)))
    spans = [20, 1, 1, 3, 6, 2, 2, 3, 1]
    # expected: [3, 1, 4, 0, 2, 0, 0, 4, 0]
    trace(0, '{}'.format(in_range(v_vect,ranges,spans)))
    # smoke and time test get_variant_array() and get_dna_ids()
    trace(0, 'get_dna_ids at {:.2f} seconds...'.format(time.time() - t0))
    ids = get_dna_ids(db)
    trace(0, 'get_variant_array at {:.2f} seconds...'.format(time.time() - t0))
    arr = get_variant_array(db, ids)
    trace(0, '...done at {:.2f} seconds'.format(time.time() - t0))
    # smoke test get_variant_defs
    trace(0, '{}'.format(get_variant_defs(db, [1,2,20,21])))
