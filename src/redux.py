#!/usr/bin/env python3
#
# Copyright (c) 2018 the Authors
#
# Purpose: Reduction and comparison main driver for Y-chromosome NGS test data
#
# Usage:
#   run script as a command with various args
#
# Environment:
#   REDUX_PATH must be set to source directory
#   config.yaml is read for configuration settings
#
import sys
import argparse
import yaml
import os
import time
from lib import *
from db import DB
from array_api import *
from sort import *


# configure default logging of any errors that occur during setup
trace = Trace(1)

# environment variable required
try:
    sys.path.append(os.environ['REDUX_PATH'])
    REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
except:
    trace(0,"Missing environment variable REDUX_PATH. Aborting.")
    sys.exit()

# parse the remainder of the configuration settings
config = yaml.load(open(REDUX_CONF))

# set up logging for diagnostics and status messages
trace = Trace(config['verbosity'])

compute_t0 = time.clock()
wall_t0 = time.time()
trace (1, 'Beginning run at [{}]'.format(time.ctime(wall_t0)))

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
parser.add_argument('-k', '--kits', help='update list of the kits from kits.txt', action='store_true')
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
parser.add_argument('-mi', '--matrix_imperf', help='matrix_imperf', action='store_true')
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

parser.add_argument('-van', '--variant_calls_name', help='variant calls name')
parser.add_argument('-vai', '--variant_calls_id', help='variant calls id')
parser.add_argument('-vap', '--variant_calls_pos', help='variant calls pos')
parser.add_argument('-vax', '--variant_calls_vix', help='variant calls vix')

parser.add_argument('-cp', '--clade_priority', help='clade priority')

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
    db.commit()

# load kits that were found in H-R web API and in zipdirs
if args.loadkits:
    db = DB(drop=False, fastload=True)
    populate_from_dataset(db)
    db.commit()

# load kits that were found in H-R web API and in zipdirs
if args.kits:
    db = DB(drop=False)
    populate_fileinfo(db, fromweb=config['use_web_api'])
    populate_analysis_kits(db)
    populate_excludes(db)
    db.commit()

# run unit tests - this is for development, test and prototyping
# not part of the actual program
if args.testdrive:
    db = db_creation()
    populate_from_dataset(db)
    trace(1,'get analysis DNA ids')
    ids = get_analysis_ids(db) # only get the kits of interest
    if len(ids) < 26:
        trace(1, 'calculate csv at {:.2f} s'.format(time.time() - wall_t0))
        out = get_variant_csv(db,ids)
        trace(1, 'write csv at {:.2f} s'.format(time.time() - wall_t0))
        open('csv.out','w').write(out)
    else:
        trace(0, 'skipping large csv file')
    # no other work flow
    trace(5, 'commit work at {:.2f} s'.format(time.time() - wall_t0))
    db.commit()
    trace(1, 'finished at {:.2f} s'.format(time.time() - wall_t0))
    trace(1, 'CPU time: {:.2f} s'.format(time.clock()))
    sys.exit(0)


# (beg) sort prototype stuff

def init_vt(initSort=False):
    vt = Variant()
    if initSort:
        vt.sort = Sort()
    return vt
    
if args.sort:
    sort = Sort()
    sort.init_sort_matrix()

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

if args.matrix_imperf:
    vt = init_vt(initSort=True)
    vt.matrix(args.rest,imperfOnly=True)

if args.unknowns:
    sort = Sort()
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
    vt = init_vt(initSort=True)
    vt.ref_name([args.variant_clade_name]+args.rest,clade=True)

if args.variant_clade_id:
    vt = init_vt(initSort=True)
    vt.ref_id([args.variant_clade_id]+args.rest,clade=True)

if args.variant_clade_pos:
    vt = init_vt(initSort=True)
    vt.ref_pos([args.variant_clade_pos]+args.rest,clade=True)

if args.variant_clade_vix:
    vt = init_vt(initSort=True)
    vt.ref_vix([args.variant_clade_vix]+args.rest,clade=True)

if args.variant_calls_name:
    vt = init_vt(initSort=True)
    vt.ref_name([args.variant_calls_name]+args.rest,calls=True)

if args.variant_calls_id:
    vt = init_vt(initSort=True)
    vt.ref_id([args.variant_calls_id]+args.rest,calls=True)

if args.variant_calls_pos:
    vt = init_vt(initSort=True)
    vt.ref_pos([args.variant_calls_pos]+args.rest,calls=True)

if args.variant_calls_vix:
    vt = init_vt(initSort=True)
    vt.ref_vix([args.variant_calls_vix]+args.rest,calls=True)

if args.variant_upd_unk:
    vt = init_vt(initSort=True)
    vt.upd_unk(args.variant_upd_unk)

if args.clade_priority:
    vt = init_vt(initSort=True)
    vt.clade_priority([args.clade_priority]+args.rest)

# (end) sort prototype stuff


trace(1, 'done at {:.2f} seconds'.format(time.time() - wall_t0))

