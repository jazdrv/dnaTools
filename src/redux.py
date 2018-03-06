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
parser.add_argument('-vu', '--variant_upd_unk', help='variant upd unknowns')
parser.add_argument('-vss', '--variant_stash', help='variant stash')
parser.add_argument('-vsu', '--variant_unstash', help='variant unstash')
parser.add_argument('-m', '--matrix', help='matrix', action='store_true')
parser.add_argument('-mp', '--matrix_perf', help='matrix_perf', action='store_true')
parser.add_argument('-u', '--unknowns', help='unknowns', action='store_true')
parser.add_argument('-o', '--sort', help='sort matrix data prototype (s_ schema currently)', action='store_true')

parser.add_argument('-vln', '--variant_lib_name', help='variant lib name reference')
parser.add_argument('-vli', '--variant_lib_id', help='variant lib id reference')
parser.add_argument('-vlp', '--variant_lib_pos', help='variant lib pos reference')
parser.add_argument('-vlx', '--variant_lib_vix', help='variant lib vix reference')

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

if args.sort:
    sort = Sort()
    sort.dbo = DB1()
    sort.sort_schema()
    sort.sort_matrix()

if args.variant_info:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.info(args.variant_info)

if args.variant_proc:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.proc(args.variant_proc)

if args.variant_proc_nonsplits:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.proc_nonsplits()

if args.variant_stash:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.stash(args.variant_stash)

if args.variant_unstash:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.unstash(args.variant_unstash)

if args.matrix:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.matrix(args.rest)

if args.matrix_perf:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.matrix(args.rest,perfOnly=True)

if args.unknowns:
    sort = Sort()
    sort.dbo = DB1()
    sort.stdout_unknowns()

if args.variant_lib_name:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.lib_name([args.variant_lib_name]+args.rest)

if args.variant_lib_id:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.lib_id([args.variant_lib_id]+args.rest)

if args.variant_lib_pos:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.lib_pos([args.variant_lib_pos]+args.rest)

if args.variant_lib_vix:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.lib_vix([args.variant_lib_vix]+args.rest)

if args.variant_upd_unk:
    vt = Variant()
    vt.dbo = DB1()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.upd_unk([args.variant_upd_unk]+args.rest)

# (end) sort prototype stuff

trace(0, "** script complete.\n")
trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))

