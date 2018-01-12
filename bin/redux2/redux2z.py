#!/usr/bin/env python3

# authors/licensing {{{

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,argparse,yaml,os,glob,shutil,re,time,csv,zipfile
from collections import defaultdict
from misc import *
from lib import *
from db import *

# }}}
# setup notes {{{

# beyond the modules specified in the lib areas, two env variables are
# required (ie: can be set in ~/.bashrc):

# REDUX_PATH - where these bin and conf files are
# REDUX_ENV - where the data will be

# }}}

# trace

start_time = time.clock()
trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))

# env

try:
    sys.path.append(os.environ['REDUX_PATH'])
    REDUX_CONF = os.environ['REDUX_PATH']+'/config.yaml'
except:
    print "Missing environment variable REDUX_CONF. Aborting."
    sys.exit()
try:
    REDUX_ENV = os.environ['REDUX_ENV']
except:
    print "Missing environment variable REDUX_ENV. Aborting."
    sys.exit()

config = yaml.load(open(REDUX_CONF))
#print ('\n'+sys.argv[0]+' version: '+config['VERSION']+'\n')

# arg parser (we can replace this -if useful- with getOpts or whatever)

parser = argparse.ArgumentParser()
parser.add_argument('-A', '--all', help='perform all possible steps', action='store_true')
parser.add_argument('-b', '--backup', help='do a backup', action='store_true')
parser.add_argument('-p', '--prep', help='prep file structure', action='store_true')
parser.add_argument('-d', '--data', help='SNP data processing', action='store_true')
args = parser.parse_args()

print ""

if args.all:
    go_backup()
    go_prep()
    go_db()
else:
    if args.backup:
        go_backup()
    if args.prep:
        go_prep()
    if args.data:
        go_db()

print "** script complete.\n"

