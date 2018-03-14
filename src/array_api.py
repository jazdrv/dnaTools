#!/usr/bin/env python3
# coding: utf-8

# Contributors: Jef Treece, Harald Alvestrand, Zak Jones, Iain McDonald
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

from collections import defaultdict
import sys, yaml

config = yaml.load(open('config.yaml'))

# diagnostics
def trace (level, msg, stream=sys.stderr):
    if level <= config['verbosity']:
        if level == 0:
            print(msg)
        else:
            print(msg, file=stream)
            stream.flush()

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


# efficiently determine if a set of positions is contained in a set of ranges
# return vector of True values if v contained in a range, else False
# does not consider the endpoints
# ranges must be sorted on their minaddr
# input vector must be sorted
def in_range(v_vect, ranges, spans):
    #v_vect = sorted(v_vect)
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


# build a dictionary of people-variants
# ppl is a 1-d array of dataset IDs we're interested in
# all entries in vcfcalls are returned; not limited to repeated calls
def get_variant_array(db, ppl):
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

# create a csv file of 2d person x variant array
def get_variant_csv(db, ppl):
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

# get the list of populated (ones that have vcfcalls entries) DNAIDs
def get_dna_ids(db):
    return [x[0] for x in db.dc.execute('select distinct(pID) from vcfcalls')]


# get build identifier by its name; creates new entry if needed
# known aliases are reduced to one entry
def get_build_byname(db, buildname='hg38'):
    if buildname.lower().strip() in ('hg19', 'grch37', 'b19', 'b37'):
        buildname = 'hg19'
    elif buildname.lower().strip() in ('hg38', 'grch38', 'b38'):
        buildname = 'hg38'
    dc = db.dc.execute('select id from build where buildNm=?', (buildname,))
    bid = None
    for bid, in dc:
        continue
    if not bid:
        dc.execute('insert into build(buildNm) values (?)', (buildname,))
        bid = dc.lastrowid
    return bid

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
    xl = list([(v[0],v[1]) for v in calls])
    iv = [v[0] for v in xl]
    pv = [v[1] for v in xl]
    trace(500, '{} calls: {}...'.format(len(pv), pv[:20]))
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
    coverage = in_range(pv, rv, spans)
    dc.close()
    rc.close()
    return iv, coverage


# find coverage for a list of kits and a list of variants
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
        trace(2, 'indels for kit {}...'.format(pid))
        trace(3, 'get_call_coverage(db, {}, {}, {})'.format(pid,[(i[0],i[1]) for i in enumerate(indel_ids)][:50],[(i[0],i[1]) for i in enumerate(spans)][:50]))
        iv,cv = get_call_coverage(db, pid, indel_ids, spans)
        trace(4, 'indels:{}..., coverage:{}...'.format([(i[0],i[1]) for i in enumerate(iv)][:50], [(i[0],i[1]) for i in enumerate(cv)][:50]))
        # store "not-covered" since it's sparse
        for cov,vid in zip(cv, iv):
            if not cov:
                cdict[pid][vid] = False
        # get snp coverage for a kit
        trace(2, 'snps for kit {}...'.format(pid))
        trace(3, 'get_call_coverage(db, {}, {})'.format(pid,snp_ids))
        iv,cv = get_call_coverage(db, pid, snp_ids)
        trace(4, 'snps:{}..., coverage:{}...'.format(iv[:5], cv[:5]))
        # store "not-covered" since it's sparse
        for cov,vid in zip(cv, iv):
            if not cov:
                cdict[pid][vid] = False

    return cdict
