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

start_time = time.clock()
t0 = time.time()
trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))

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
    trace(0,'get analysis DNA ids')
    ids = get_analysis_ids(db) # only get the kits of interest
    if len(ids) < 26:
        trace(0, 'calculate csv at {}'.format(time.clock()))
        out = get_variant_csv(db,ids)
        trace(0, 'write csv at {}'.format(time.clock()))
        open('csv.out','w').write(out)
    else:
        trace(0, 'skipping large csv file')
    # no other work flow
    trace(0, 'commit work at {}'.format(time.clock()))
    db.commit()
    trace(1, 'commit done at {:.2f} seconds'.format(time.time() - t0))
    db.close()
    sys.exit(0)


trace(0, "** script complete.\n")
trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))

