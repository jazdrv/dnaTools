#!/usr/bin/env python3

# Copyright (c) 2018 the Authors

# Contributors: Jef Treece, Harald Alvestrand, Iain McDonald, Zak Jones
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html


import sys
import argparse
import yaml
import os
import glob
import shutil
import re,time,csv,zipfile
from collections import defaultdict
from lib import *
from db import DB
from array_api import *

# (beg) sort prototype libs

from db1 import DB1 #this one can be better merged with db.py (but for debugging, I prefer this one) 
from sort import *

# (end) sort prototype libs

# required environment vars:
# REDUX_PATH - where the source code and config.yaml lives

start_time = time.clock()
t0 = time.time()
trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))

# environment variable required
#try:
#    sys.path.append(os.environ['REDUX_PATH'])
#    REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
#except:
#    trace(0,"Missing environment variable REDUX_PATH. Aborting.")
#    sys.exit()

# parse the remainder of the configuration settings
#config = yaml.load(open(REDUX_CONF))

try:
    config = yaml.load(open(os.environ['REDUX_CONF_ZAK']))
except:
    print("Missing environment variable REDUX_CONF_ZAK. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

# basic strategy for command-line arguments
#  -command-line args mainly select one or more basic execution elements (below)
#  -these are standalone work units
#  -config settings (above) are for almost all of the fine tuning
#  -limited number of options from the command line

# basic execution elements selected by command line arguments
#  -running: logging, help messages
#  -tasks: workflows, such as compute a particular type of output
#  -maintenance: backup, cleanup, diagnostics, tools
#  -output: controlling type of output files


# arg parser
parser = argparse.ArgumentParser()
# running

# tasks
parser.add_argument('-a', '--all', help='perform all possible steps (prob best not to use for now)', action='store_true')
parser.add_argument('-c', '--create', help='clean start with a new database', action='store_true')
parser.add_argument('-l', '--loadkits', help='load all of the kits', action='store_true')
parser.add_argument('-t', '--testdrive', help='runs some unit tests', action='store_true')

# (beg) sort prototype stuff

parser.add_argument('-vi', '--variant_info', help='variant info', type=str)
parser.add_argument('-vp', '--variant_proc', help='variant process', type=str)
parser.add_argument('-vpns', '--variant_proc_nonsplits', help='variant process all non splits', action='store_true')
parser.add_argument('-vu', '--variant_upd_unk', help='variant upd unknowns', type=str)
parser.add_argument('-vss', '--variant_stash', help='variant stash')
parser.add_argument('-vsu', '--variant_unstash', help='variant unstash')
parser.add_argument('-m', '--matrix', help='matrix', action='store_true')
parser.add_argument('-mp', '--matrix_perf', help='matrix_perf', action='store_true')
parser.add_argument('-u', '--unknowns', help='unknowns', action='store_true')
parser.add_argument('-o', '--sort', help='sort matrix data prototype (s_ schema currently)', action='store_true')

parser.add_argument('-vrn', '--variant_ref_name', help='variant ref name')
parser.add_argument('-vri', '--variant_ref_id', help='variant ref id')
parser.add_argument('-vrp', '--variant_ref_pos', help='variant ref pos')
parser.add_argument('-vrx', '--variant_ref_vix', help='variant ref vix')

parser.add_argument('-vcn', '--variant_clade_name', help='variant clade name')
parser.add_argument('-vci', '--variant_clade_id', help='variant clade id')
parser.add_argument('-vcp', '--variant_clade_pos', help='variant clade pos')
parser.add_argument('-vcx', '--variant_clade_vix', help='variant clade vix')

parser.add_argument('rest', nargs=argparse.REMAINDER)

# (end) sort prototype stuff

# maintenance
parser.add_argument('-b', '--backup', help='do a "backup"', action='store_true')

# output

args = parser.parse_args()




# main program

# drop database and have a clean start
if args.create:
    db = db_creation()

# load kits that were found in H-R web API and in zipdirs
if args.loadkits:
    populate_from_dataset(db)

# run unit tests - this is for development, test and prototyping
# not part of the actual program
if args.testdrive:
    db = db_creation()
    populate_from_dataset(db)
    trace(0,'get DNA ids')
    ids = get_dna_ids(db)
    if config['kitlimit'] < 26:
        trace(0, 'calculate array')
        out = get_variant_csv(db,ids)
        open('csv.out','w').write(out)
    else:
        trace(0, 'skipping large csv file')
    # no other work flow
    trace(0, 'commit work')
    db.commit()
    db.close()
    sys.exit(0)

# run everything
if args.all:
    go_backup()
    go_prep()
    go_db()


# (beg) sort prototype stuff

def init_vt(initSort=False):
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    if initSort:
        vt.sort = Sort()
        vt.sort.dbo = vt.dbo
    return vt

def init_sort():
    sort = Sort()
    sort.dbo = DB1()
    return sort

if args.sort:
    sort = init_sort()
    sort.sort_schema()
    sort.sort_matrix()

if args.variant_info:
    vt = init_vt(initSort=True)
    vt.info(args.variant_info)

if args.variant_proc:
    vt = init_vt(initSort=True)
    vt.proc(args.variant_proc)

if args.variant_proc_nonsplits:
    vt = init_vt(initSort=True)
    vt.proc_nonsplits()

if args.variant_stash:
    vt = init_vt(initSort=True)
    vt.stash(args.variant_stash)

if args.variant_unstash:
    vt = init_vt(initSort=True)
    vt.unstash(args.variant_unstash)

if args.matrix:
    vt = init_vt(initSort=True)
    vt.matrix(args.rest)

if args.matrix_perf:
    vt = init_vt(initSort=True)
    vt.matrix(args.rest,perfOnly=True)

if args.unknowns:
    sort = init_sort()
    sort.stdout_unknowns()

if args.variant_ref_name:
    vt = init_vt()
    vt.ref_name([args.variant_ref_name]+args.rest)

if args.variant_ref_id:
    vt = init_vt()
    vt.ref_id([args.variant_ref_id]+args.rest)

if args.variant_ref_pos:
    vt = init_vt()
    vt.ref_pos([args.variant_ref_pos]+args.rest)

if args.variant_ref_vix:
    vt = init_vt()
    vt.ref_vix([args.variant_ref_vix]+args.rest)

if args.variant_clade_name:
    vt = init_vt()
    vt.ref_name([args.variant_clade_name]+args.rest,clade=True)

if args.variant_clade_id:
    vt = init_vt()
    vt.ref_id([args.variant_clade_id]+args.rest,clade=True)

if args.variant_clade_pos:
    vt = init_vt()
    vt.ref_pos([args.variant_clade_pos]+args.rest,clade=True)

if args.variant_clade_vix:
    vt = init_vt()
    vt.ref_vix([args.variant_clade_vix]+args.rest,clade=True)

if args.variant_upd_unk:
    vt = init_vt(initSort=True)
    vt.upd_unk(args.variant_upd_unk)

# (end) sort prototype stuff

trace(0, "** script complete.\n")
trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))

