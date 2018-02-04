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

# required environment vars:
# REDUX_PATH - where the source code and config.yaml lives

start_time = time.clock()
t0 = time.time()
trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))

# environment variable required
try:
    sys.path.append(os.environ['REDUX_PATH'])
    REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
except:
    trace(0,"Missing environment variable REDUX_PATH. Aborting.")
    sys.exit()

# parse the remainder of the configuration settings
config = yaml.load(open(REDUX_CONF))


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
    if config['kitlimit'] < 15:
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


trace(0, "** script complete.\n")
trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))

