#!/usr/bin/env python3
#
# Copyright (c) 2018 the Authors
#
# Purpose: Reduction and comparison array and sort
#
# Usage:
#   import Variant or Sort class or run directly as a script for unit tests
#
import sys,os,yaml,csv,json,time,numpy as np
from beautifultable import BeautifulTable
from collections import OrderedDict
import pandas as pd
from array_api import *
import pickle
from db import DB
from lib import Trace
import shelve

REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
config = yaml.load(open(REDUX_CONF))

trace = Trace(level=config['verbosity'])


def l2s(lst):
    return ",".join(str(x) for x in lst)

# Class: VKcalls
# Purpose: a data structure to store an array of calls and, helper methods
# Input:
#   calls_by_vid: associative array, key is a vid, val is a vector of calls
#   kvect: pid corresponding to enries in vector of calls (kits)
class VKcalls(object):
    def __init__(self, dbo, calls_by_vid, kvect):
        self.dbo = dbo
        # vector of variants (vIDs)
        self.vids = list(calls_by_vid.keys())
        # number of variants
        self.nvids = len(self.vids)
        # variant definitions
        self.vdefs = None
        # vector of kits (pIDs) corresponding to calls per vID
        self.kits = kvect
        # number of kits
        self.nkits = len(kvect)
        # kit definitions (so we can display kitID)
        self.kdefs = None
        # kit order - vector of indices into self.kits
        self.korder = range(self.nkits)
        # variant order - vector of indices into self.vids
        self.vorder = range(self.nvids)
        # call information (0=ref, 1=alt, -1=nc, -2=mixed)
        # pandas matrix in (row,col) = (variant,kit) order
        self.calls = pd.DataFrame(calls_by_vid).T.rename(
                         dict(zip(range(self.nkits),self.kits)), axis=1)
        trace(5, 'calls:\n{}'.format(self.calls))

    # Method: update_data
    # Purpose: replace data with a new array (e.g. a sorted one)
    # Info:
    #   does not update vdefs,kdefs, which is an issue if kits changed
    def update_data(self, A):
        self.calls = A
        self.kits = list(A.axes[1])
        self.nkits = len(self.kits)
        self.korder = list(range(self.nkits))
        self.vids = list(A.axes[0])
        self.nvids = len(self.vids)
        self.vorder = list(range(self.nvids))

    # represent this data structure as a human-readable string
    def __str__(self):
        if not self.kdefs:
            self.kdefs = get_kit_ids(self.dbo, [self.kits[ii] for ii in
                                                    self.korder])
        NP = self.calls.as_matrix()
        table = BeautifulTable(max_width=450)
        table.column_headers = ['vID']+[str(self.kdefs[self.kits[ii]])
                                            for ii in self.korder]
        table.append_row(['pID-->']+[self.kits[ii] for ii in self.korder])
        table.append_row(['']+['']*self.nkits)
        for iV in self.vorder:
            V = self.vids[iV]
            rowvals = NP[iV].tolist()
            table.append_row([V]+[rowvals[ii] for ii in self.korder])
            table.row_seperator_char = ''
            table.column_seperator_char = ''
            table.column_alignments['vID'] = BeautifulTable.ALIGN_LEFT
        trace(1,'Total matrix kits: {}'.format(self.nkits))

        return '{}'.format(table)

    # represent this data structure as a csv
    def to_csv(self, fname):
        if not self.kdefs:
            self.kdefs = get_kit_ids(self.dbo, [self.kits[ii] for ii in
                                                    self.korder])
        if not self.vdefs:
            vd = get_variant_defs(self.dbo, self.vids)
            self.vdefs = dict([(v[0],(v[1],v[2],v[3])) for v in vd])
        NP = self.calls.as_matrix()
        fieldnames = ['vID','pos','anc','der'] +\
                [str(self.kdefs[self.kits[ii]]) for ii in self.korder]
        with open (fname, 'w') as csvfile:
            cf = csv.DictWriter(csvfile, fieldnames=fieldnames)
            cf.writeheader()
            d = ['','','',''] + [self.kits[ii] for ii in self.korder]
            d = dict(zip(fieldnames,d))
            cf.writerow(d)
            for iV in self.vorder:
                V = self.vids[iV]
                trace(2,'{}'.format(V))
                pos,ref,alt = self.vdefs[V]
                rowvals = [V, pos, ref, alt]
                cvals =  NP[iV].tolist()
                rowvals += [cvals[ii] for ii in self.korder]
                rowdict = dict(zip(fieldnames,rowvals))
                trace(2,'{}'.format(rowdict))
                cf.writerow(rowdict)
    
    # Method: partition4
    # Purpose: partition array into four quadrants
    #   A - upper left, contains the largest clade
    #   B - upper right, complement set for A (zeroes)
    #   C - lower left, variants not in largest clade
    #   D - clades parallel to A
    #   E - any variants where all kits are called
    # Info:
    #   sort on descending call count by variant - highest number defines A;
    #   kits not in A should have all zeroes for these variants -> kits define
    #   B variants not in A should have all zeroes for kits in A -> defines C
    #   variants not in A should have non-zero for some kits in B -> defines D
    #   in case there's a variant called by all kits, B and C may be empty
    #
    #       c1       c2
    #    11111111|0000000
    #      ...   |   ...  
    #r1     A    |    B   
    #     -------+--------
    #    00000000|1100101
    #      ...   |   ... 
    #r2     C    |    D
    #    
    def partition4(self, A):

        # count of kits that have a call for the variant
        sums = A.sum(axis=1)
        idxsums = list(zip(sums.T.axes[0], sums.values))

        # if there are variants where all kits are called, they go in E
        fullrows = [x[0] for x in idxsums if x[1] == len(A.columns)]
        nonfullrows = [x[0] for x in idxsums if x[1] != len(A.columns)]

        E = A.filter(fullrows, axis=0)

        # remove variants (rows) from A if we put them in E
        A = A.filter(nonfullrows, axis=0)

        # partition the set of kits on the variant with the most calls
        sums = A.sum(axis=1)
        idx = sums.idxmax()
        colidx = list(A.T.axes[0])
        colfilt = list(A.T[idx])
        # keep indices where call exists (value > 1) --> A and C together
        grp1kits = [x[1] for x in zip(colfilt, colidx) if x[0] > 0]
        grp2kits = [x[1] for x in zip(colfilt, colidx) if x[0] == 0]
        AC = A.filter(grp1kits, axis=1)

        # partition on variants that have kits under these calls (sum>0)
        sums = AC.sum(axis=1)
        idxsums = list(zip(sums.T.axes[0], sums.values))
        grp1rows = [x[0] for x in idxsums if x[1] > 0]
        grp2rows = [x[0] for x in idxsums if x[1] == 0]

        return fullrows, grp1rows, grp2rows, grp1kits, grp2kits

    # Method: partition_sort
    # Purpose: partition an array into clades
    # Recursive algorithm:
    #   Input array A
    #   partition A into five sub-arrays
    #     E=snps where all kits have the call
    #     A=largest clade (highest count of calls across kits)
    #     B=parallel clade to A
    #     C=complement kits to those to D
    #     D=minor clade or set of clades
    #   call partition_sort on A and D
    #   reindex A with the new values from sorting A and D
    # Info:
    #   this should work well for perfect variants (those that are either
    #   positively or negatively called for all kits in the set, with no kits
    #   ambiguous) and it may fail if there are recurrent variants and/or
    #   ambiguous calls
    # WARN: recursive, may not scale
    def partition_sort(self, A):

        # don't recurse forever - this may need some work
        if A.shape[0] < 2 or A.shape[1] < 2:
            return A

        rE,r1,r2,c1,c2 = self.partition4(A)
        Ap = A.filter(r1, axis=0).filter(c1, axis=1)
        Dp = A.filter(r2, axis=0).filter(c2, axis=1)

        Anew = self.partition_sort(Ap)
        r1A, c1A = list(Anew.axes[0]), list(Anew.axes[1])
        Dnew = self.partition_sort(Dp)
        r2D, c2D = list(Dnew.axes[0]), list(Dnew.axes[1])

        rows = rE+r1A+r2D
        cols = c1A+c2D

        trace(10, 'cols: {}'.format(c1A+c2D))
        return A.reindex(rows).reindex(columns=cols)


# Class: Sort
# Purpose: container and methods for sorting variants and kits
class Sort(object):
    def __init__(self):
        self.dbo = DB(drop=False)
        self.KITS = None
        self.VARIANTS = None
        self.NPu = None # unknowns matrix
        self.NPp = None # perfect calls matrix
        self.NPi = None # identical calls matrix
        self.NPm = None # mixed calls matrix
        self.SS = None #sups+subs
        self.identcallDATA = {}
        self.perfectDATA = {}
        self.unknownDATA = {}
        self.mixedDATA = {}

    def sort_matrix(self):
        #get data
        self.create_mx_data()

        #sort

        #cache sups + subs
        # self.reset_ss_data()

        #save data
        self.save_mx()

    def sort_schema(self):
        #Note: where the sort DDL schema is first run
        self.dbo.create_schema(schemafile='sort-schema.sql')

    def save_mx(self):
        trace(1, 'beg MatrixData stash at {}'.format(time.clock()))
        with shelve.open('saved-data') as db:
            db['KITS'] = self.KITS
            db['VARIANTS'] = self.VARIANTS
            db['perfectDATA'] = self.perfectDATA
            db['mixedDATA'] = self.mixedDATA
            db['identcallDATA'] = self.identcallDATA
            db['unknownDATA'] = self.unknownDATA
            db['NPp'] = self.NPp
            db['vdefs'] = self.vdefs

    def restore_mx_data(self):
        trace(1, 'beg MatrixData restore at {}'.format(time.clock()))
        with shelve.open('saved-data') as db:
            self.KITS = db['KITS']
            self.VARIANTS = db['VARIANTS']
            self.perfectDATA = db['perfectDATA']
            self.mixedDATA = db['mixedDATA']
            self.identcallDATA = db['identcallDATA']
            self.unknownDATA = db['unknownDATA']
            self.NPp = db['NPp']
            self.vdefs = db['vdefs']

    # Method: create_mx_data
    # Purpose: create an array of genotypes out of call data
    # Side effect: matrices are populated in class variables
    # Info:
    #   matrix[vid] is a vector of genotype IDs for the given variant
    #     the index into this vector corresponds to kitid
    #     a genotype of -1 means unknown due to poor coverage
    #     a genotype of -2 means gt was ambiguous (0/1, 0/2, 1/2, etc)
    def create_mx_data(self):

        print("beg MatrixData create: %s" % format(time.clock()))

        # get all call info (arr) and coverage info (cov)
        ppl = get_analysis_ids(self.dbo)
        arr, ppl, vids = get_variant_array(self.dbo, ppl)
        trace(5, 'arr: {}...'.format(arr[ppl[0]]))
        cov = get_kit_coverages(self.dbo, ppl, vids)

        # loop through (vid,refid,altid) for our vids of interest
        # value in resulting matrix:
        #   = allele ID if genotype is solidly called
        #   = -1 if position is not covered
        #   = -2 if the call is mixed (ref or alt possible)
        # the ancestral value is always the allele id of the variant anc
        # values are kept as database IDs for efficiency
        # can be converted to human-readable later when display is needed

        # called variants that are identical across all kits
        identcallDATA = {}
        # variants that are called for every kit and not all kits are the same
        perfectDATA = {}
        # variants that have no solid calls
        unknownDATA = {}
        # variants that have some solid calls and some unknowns
        mixedDATA = {}

        vd = get_variant_defs(self.dbo, vids)
        trace(2, 'vdefs({}): [{},{}....{},{}]'.format(len(vd),
                                    vd[0],vd[1], vd[-2],vd[-1]))
        self.vdefs = dict([(v[0],(v[1],v[2],v[3])) for v in vd])

        # main loop to populate the matrices
        for vv in self.vdefs:
            vect = []
            for pp in ppl:
                gt = None
                try:
                    # there is a call for this pid,vid
                    pf,igt = arr[pp][vv]
                    covered = pf
                    if igt in (0,1):
                        # gt from get_variant_array is 0,1 ==> anc,der
                        gt = igt
                    else:
                        # gt is not 0,1 ==> ambiguous call: don't know
                        gt = -2
                except:
                    # no call for this pid,vid
                    try:
                        # if variant is not covered, we don't know genotype
                        if cov[pp][vv] == 0:
                            covered = False
                            gt = -1
                        elif cov[pp][vv] in (1,2,3):
                            # covered, gt=anc by get_variant_array convention
                            covered = True
                            gt = 0
                        else:
                            trace(0, 'unexpected value for coverage')
                    # no entry in coverage is default (covered)
                    except:
                        covered = True
                        gt = 0
                if gt is None:
                    raise ValueError('gt not set')
                vect.append(gt)
            trace(3,'vect: {}...'.format(vect[:10]))

            # no unknowns
            if -1 not in vect and -2 not in vect:
                # we have at least one ref and one alt
                if (0 in vect) and (1 in vect):
                    perfectDATA[vv] = vect
                # all the same genotype across kits
                elif vect.count(vect[0]) == len(vect):
                    identcallDATA[vv] = vect
                else:
                    raise ValueError('unhandled vector type')
            # some unknowns exist (due to no coverage or mixed calls)
            else:
                # there's at least one with a call
                if (len([t for t in vect if t>=0]) > 0):
                    mixedDATA[vv] = vect
                else:
                    unknownDATA[vv] = vect

        # FIXME - used downstream, this definition isn't correct yet
        self.VARIANTS = vids
        self.KITS = ppl
        trace(2,'KITS: {}'.format(self.KITS))
        trace(2,'VARIANTS: {}...'.format(self.VARIANTS[:20]))

        # show some debugging info
        self.NPp = np.matrix(list(perfectDATA.values()))
        trace(2,'NPp:\n{}'.format(self.NPp))
        self.NPi = np.matrix(list(identcallDATA.values()))
        trace(2,'NPi:\n{}'.format(self.NPi))
        self.NPm = np.matrix(list(mixedDATA.values()))
        trace(2,'NPm:\n{}'.format(self.NPm))
        self.NPu = np.matrix(list(unknownDATA.values()))
        trace(2,'NPu:\n{}'.format(self.NPu))

        # summarize the results for diagnostics
        trace(1,'Total matrix kits: {}'.format(len(self.KITS)))
        trace(1,'Total perfect matrix variants: {}'.format(len(self.NPp)))
        trace(1,'Total identical matrix variants: {}'.format(len(self.NPi)))
        trace(1,'Total mixed matrix variants: {}'.format(len(self.NPm)))
        trace(1,'Total unknowns matrix variants: {}'.format(len(self.NPu)))
        trace(1,'Total matrix variants: {}'.format(len(self.VARIANTS)))

        self.identcallDATA = identcallDATA
        self.perfectDATA = perfectDATA
        self.unknownDATA = unknownDATA
        self.mixedDATA = mixedDATA
        self.save_mx()


# Class: Variant
# Purpose: not sure, seems like an empty shell for Sort class
class Variant(Sort):

    def __init__(self):
        Sort.__init__(self)

    def proc(self,vname):
        #Note: just process one variant 
        self.restore_mx_data()
        self.vix = self.proc_vname(vname)
        if self.vix is None:
            print ("No variants could be found. Exiting.")
            sys.exit()
        self.proc_chk(allowImperfect=False)
        
    def info(self,vname):
        self.restore_mx_data()
        allowImperfect = config['allowImperfectWithVInfo']
        self.vix = self.proc_vname(vname)
        self.set_info(lev=2,allowImperfect=allowImperfect)
        self.stdout_info()
        
    def matrix(self,argL=None):

        self.restore_mx_data()
        kix = self.KITS
        #m = VKcalls(self.dbo, self.mixedDATA, kix)
        #trace(1,'{}'.format(m))
        #m = VKcalls(self.dbo, self.unknownDATA, kix)
        #trace(1,'{}'.format(m))
        #m = VKcalls(self.dbo, self.identcallDATA, kix)
        #trace(1,'{}'.format(m))
        m = VKcalls(self.dbo, self.perfectDATA, kix)
        A = m.partition_sort(m.calls)
        m.update_data(A)
        m.to_csv('csv.out')
        trace(3,'m:\n{}'.format(m))


# test framework
if __name__=='__main__':
    v = Variant()
    try:
        v.matrix()
    except:
        raise
        v.create_mx_data()
        v.matrix()
        v.to_csv('csv.out')
