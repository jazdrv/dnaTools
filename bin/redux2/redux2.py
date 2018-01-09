#!/bin/python

# legacy redux.bash steps {{{

# -b [backup]

# (1) BACKUPS

# -p [prep]

# (2) CHECK NECESSARY FILES EXIST
# (3) UNZIP FILES
# (4) EXTRACT RAW DATA

# -d [data]

# (5) POST-PROCESSING OF RAW RESULTS FOR ODDITIES
# (6) DATA CATEGORISATION (FORM TABLES)
# (7) DATA ORDERING
# (8) NAMES TO SNPS
# (9) STATISTICS GENERATION

# [form+tree]

# (10) FORM REPORT
# (11) TREE STRUCTURE

# -q [qc] 

# (12) CONSISTENCY CHECKING
# (13) FLAG CLADES THAT CAN BE MERGED
# (14) IDENTIFY FORCED POSITIVES : CLADES THAT SHOULDN'T BE MERGED?

# -r [rpts]

# (15) INSERT TREE INTO REPORT
# (16) SHORT REPORT
# (17) HTML REPORT

# -a [ages/tmrca]

# (18) INITIAL AGE ANALYSIS

# -o [output]

# (19) TREE MERGING
# (20) AGE REPORT WRITING
# (21) TREE SVG WRITING

# }}}

import sys,argparse, yaml

def go_all():
    go_backup()
    go_prep()
    go_data()
    
def go_backup():
    print "** performing backup."
    print "** + backup done."
    
def go_prep():
    print "** prepare file structure."
    print "** + prep done."
    
def go_data():
    print "** process SNP data."
    print "** + SNP processing done."

config=yaml.load(open('config.yaml'))

parser = argparse.ArgumentParser()
parser.add_argument('-A', '--all', help='perform all possible steps', action='store_true')
parser.add_argument('-b', '--backup', help='do a backup', action='store_true')
parser.add_argument('-p', '--prep', help='prep file structure', action='store_true')
parser.add_argument('-d', '--data', help='SNP data processing', action='store_true')
args = parser.parse_args()

#print ('\n'+sys.argv[0]+' version: '+config['VERSION']+'\n')
print ""

if args.all:
    go_all()
else:
    if args.backup:
        go_backup()
    if args.prep:
        go_prep()
    if args.data:
        go_data()

print "** script complete.\n"

