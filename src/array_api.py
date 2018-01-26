#!/usr/bin/env python3
# coding: utf-8

# Contributors: Jef Treece, Harald Alvestrand, Zak Jones, Iain McDonald
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html


# build a dictionary of people-variants
# ppl is a 1-d array of dataset IDs we're interested in
def variant_array(db, ppl):
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
def variant_csv(db, ppl):
    arr, pl, vl = variant_array(db, ppl)
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
