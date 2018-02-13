#!/usr/bin/env python3

# license {{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

import sys,argparse,yaml,os
from sort import *
#from db import *
from db_s import *

# }}}

try:
    config = yaml.load(open(os.environ['REDUX_CONF']))
except:
    print("Missing environment variable REDUX_CONF. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

#parser

parser = argparse.ArgumentParser()
parser.add_argument('-vi', '--variant_info', help='variant info', type=str)
parser.add_argument('-vp', '--variant_proc', help='variant process', type=str)
parser.add_argument('-m', '--matrix', help='matrix', action='store_true')
parser.add_argument('-u', '--unknowns', help='unknowns', action='store_true')
parser.add_argument('-o', '--sort', help='sort matrix data prototype (s_ schema currently)', action='store_true')
args = parser.parse_args()

#arg handling

if args.sort:
    sort = Sort()
    sort.dbo = DB()
    sort.sort_schema()
    sort.sort_ins_sample_data()
    sort.sort_matrix()

if args.variant_info:
    vt = Variant()
    vt.dbo = DB()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.info(args.variant_info)

if args.variant_proc:
    vt = Variant()
    vt.dbo = DB()
    vt.dbo.db = vt.dbo.db_init()
    vt.dbo.dc = vt.dbo.cursor()
    vt.sort = Sort()
    vt.sort.dbo = vt.dbo
    vt.proc(args.variant_proc)

if args.matrix:
    sort = Sort()
    sort.dbo = DB()
    sort.stdout_matrix(refreshDbFlg=True)

if args.unknowns:
    sort = Sort()
    sort.dbo = DB()
    sort.stdout_unknowns()
