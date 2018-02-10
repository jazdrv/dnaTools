#!/usr/bin/env python3
# coding: utf-8

# Contributors: Jef Treece, Harald Alvestrand, Zak Jones, Iain McDonald
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# build a dictionary of people-variants
# ppl is a 1-d array of dataset IDs we're interested in
def get_variant_array(db, ppl):
    from collections import defaultdict
    c1 = db.dc
    c1.execute('drop table if exists tmpt')
    c1.execute('create temporary table tmpt(id integer)')
    c1.executemany('insert into tmpt values(?)', [(p,) for p in ppl])
    c1 = c1.execute('''select distinct v.id,c.pid from variants v
                     inner join vcfcalls c on c.vid=v.id
                     inner join tmpt t on c.pid=t.id''')
    # build a 2d array
    arr = {}
    var = set()
    for pp in ppl:
        arr[pp] = defaultdict()

    for (v,p) in c1:
        arr[p][v] = 1
        var.add(v)
    c1.execute('drop table tmpt')

    return arr, ppl, list(var)

# create a csv file of 2d person x variant array
def get_variant_csv(db, ppl):
    arr, pl, vl = get_variant_array(db, ppl)
    out = ',' + ','.join(['{}'.format(v) for v in ppl]) + '\n'
    for vv in vl:
        out += '{},'.format(vv)
        for pp in pl:
            try:
                out += '{},'.format(arr[pp][vv])
            except:
                out += ','
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

# get coverage for indels
# pid is the dataset id
def get_indel_coverage(db, pid):
    from lib import get_call_coverage
    # fixme - get and reuse list of indels since it takes time to find them
    # using this query and we need to call this on each kit
    c = db.cursor()
    c = c.execute('''select distinct v.id,max(length(a.allele),length(b.allele))
                     from variants v
                     inner join alleles a on a.id=v.anc
                     inner join alleles b on b.id=v.der
                     inner join vcfcalls c on c.vid=v.id
                     where length(a.allele)>1 or length(b.allele)>1''')
    indels = [t for t in c]
    ids = list([t[0] for t in indels])
    spans = list([int(t[1]) for t in indels])
    c.close()
    cv = get_call_coverage(db, pid, ids, spans)
    return cv
