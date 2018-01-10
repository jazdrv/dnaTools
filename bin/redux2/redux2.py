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

import sys,argparse,yaml,os,glob,shutil

def go_all():
    go_backup()
    go_prep()
    go_data()
    
def go_backup():

    print "** performing backup."
    
    #BACKUPFILES = "variant-*.txt snps_hg19.csv snp-names.csv snp-used.csv report.csv short-report.csv clades.csv tree.txt raw-ages.txt"
    print "** (msg) CREATING BACKUP COPIES OF EXISTING FILES..."

    # autobackup folder
    
    if (os.path.isdir("autobackup")):
        files = glob.glob('autobackup/*')
        for f in files:
            os.remove(f)
    else:
        os.makedirs("autobackup")
    
    for FILE_PATTERN in config['backup_files'].split():
        for FILE in glob.glob(FILE_PATTERN):
            shutil.copy(FILE, 'autobackup')

    # autobackup2 folder {{{

    # Make further backup copies when running the script from scratch
    # This is useful when you want to make changes to the bad/inconsistent list, but still want to compare to the original previous run.
    # For example:
    # gawk 'NR==FNR {c[$5]++;next};c[$5]==0' tree.txt autobackup2/tree.txt
    # will tell you the changes to the tree structure that have resulted from the addition of new kits between "from-scratch" runs.

    #if (os.path.isdir("autobackup2")):
    #    files = glob.glob('autobackup2/*')
    #    for f in files:
    #        os.remove(f)
    #else:
    #    os.makedirs("autobackup2")
    
    #for FILE_PATTERN in config['backup_files'].split():
    #    for FILE in glob.glob(FILE_PATTERN):
    #        shutil.copy(FILE, 'autobackup2')

    # }}}

    if config['make_report'] > 0:
        print "MAKING REPORT..."
        if os.path.exists("report.csv"):
            os.remove("report.csv")
    
    # if config['skip_zip'] == 0:
    #   ... don't know what this is for!

    print "** + backup done."

def go_prep():
    print "** prepare file structure."
    print "** + prep done."
    
def go_data():
    print "** process SNP data."
    print "** + SNP processing done."

REDUX_CONF = os.environ['REDUX_CONF']
config = yaml.load(open(REDUX_CONF))

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

