# !/bin/python

# authors/licensing{{{

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,os,sqlite3,yaml,time,csv,numpy as np
#from misc import *
#from lib import *

# }}}

REDUX_CONF = os.environ['REDUX_CONF']
REDUX_ENV = os.environ['REDUX_ENV']
config = yaml.load(open(REDUX_CONF))
start_time = 0 # need to fix this

#note: this lib could easily be a class

def db_init(trace):
    trace (1, "Initialising database...")
    db = sqlite3.connect('variant.db')
    dc = db.cursor()
    return dc
    
def db_drop_tables(dc):
    if (config['drop_tables'] > 0):
        dc.execute('''DROP TABLE IF EXISTS variants''')
        dc.execute('''DROP TABLE IF EXISTS hg19''')
        dc.execute('''DROP TABLE IF EXISTS hg38''')
        dc.execute('''DROP TABLE IF EXISTS kits''')
        dc.execute('''DROP TABLE IF EXISTS people''')
        dc.execute('''DROP TABLE IF EXISTS strs''')
        dc.execute('''DROP TABLE IF EXISTS calls''')
        dc.execute('''DROP TABLE IF EXISTS tree''')
    
def db_create_tables(dc):

    # variants tbl

    dc.execute('''CREATE TABLE IF NOT EXISTS variants
                  (id INTEGER PRIMARY KEY,
                  grch37 INTEGER,
                  length SMALLINT,
                  ref TEXT,
                  alt TEXT,
                  inage BOOLEAN)''')

    # hg19

    dc.execute('''CREATE TABLE IF NOT EXISTS hg19
                  (id INTEGER PRIMARY KEY,
                  snp BOOLEAN,
                  grch37 INTEGER,
                  grch37end INTEGER,
                  name CHARACTER(32),
                  anc CHARACTER(64),
                  der CHARACTER(64))''')

    # hg38

    dc.execute('''CREATE TABLE IF NOT EXISTS hg38
                  (id INTEGER PRIMARY KEY,
                  snp BOOLEAN,
                  grch38 INTEGER,
                  grch38end INTEGER,
                  name CHARACTER(32),
                  anc CHARACTER(64),
                  der CHARACTER(64))''')

    # kits

    dc.execute('''CREATE TABLE IF NOT EXISTS kits
                  (id INTEGER PRIMARY KEY,
                  vcffile CHARACTER(256),
                  bedfile CHARACTER(256),
                  company CHARACTER(16),
                  test CHARACTER(16),
                  kitnum CHARACTER(16),
                  testdate INTEGER,
                  uploaddate INTEGER,
                  personid INTEGER)''')

    # people

    dc.execute('''CREATE TABLE IF NOT EXISTS people
                  (person INTEGER PRIMARY KEY,
                  ngstest BOOLEAN,
                  ftdnaid CHARAACTER(16) REFERENCES kits(kitnum),
                  fgcid CHARAACTER(16) REFERENCES kits(kitnum),
                  yseqid CHARAACTER(16) REFERENCES kits(kitnum),
                  otherid CHARAACTER(16) REFERENCES kits(kitnum),
                  mdkaname CHARACTER(16),
                  mdkacountry CHARACTER(2),
                  mdkalat REAL,
                  mdkalong REAL,
                  mdkalocweight REAL,
                  mdkadesc CHARACTER(255),
                  mdkayear SMALLINT,
                  mdkaunc SMALLINT,
                  mdkaafter SMALLINT,
                  mdkabefore SMALLINT,
                  mrkhid INTEGER REFERENCES tree(id),
                  coverage INTEGER,
                  agecoverage INTEGER,
                  totcoverage INTEGER,
                  nsnps INTEGER,
                  snpsintree SMALLINT,
                  singletons SMALLINT,
                  qsingletons SMALLINT,
                  agesingletons SMALLINT)''')

    # strs

    dc.execute('''CREATE TABLE IF NOT EXISTS strs
                  (person INTEGER PRIMARY KEY,
                  DYS393 TINYINT,
                  DYS390 TINYINT,
                  DYS19 TINYINT,
                  DYS391 TINYINT,
                  DYS385a TINYINT,
                  DYS385b TINYINT,
                  DYS426 TINYINT,
                  DYS388 TINYINT,
                  DYS439 TINYINT,
                  DYS389i TINYINT,
                  DYS392 TINYINT,
                  DYS389ii TINYINT,
                  DYS458 TINYINT,
                  DYS459a TINYINT,
                  DYS459b TINYINT,
                  DYS455 TINYINT,
                  DYS454 TINYINT,
                  DYS447 TINYINT,
                  DYS437 TINYINT,
                  DYS448 TINYINT,
                  DYS449 TINYINT,
                  DYS464a TINYINT,
                  DYS464b TINYINT,
                  DYS464c TINYINT,
                  DYS464d TINYINT,
                  DYS460 TINYINT,
                  YH4 TINYINT,
                  YCAIIa TINYINT,
                  YCAIIb TINYINT,
                  DYS456 TINYINT,
                  DYS607 TINYINT,
                  DYS576 TINYINT,
                  DYS570 TINYINT,
                  CDYa TINYINT,
                  CDYb TINYINT,
                  DYS442 TINYINT,
                  DYS438 TINYINT,
                  DYS531 TINYINT,
                  DYS578 TINYINT,
                  DYF395S1a TINYINT,
                  DYF395S1b TINYINT,
                  DYS590 TINYINT,
                  DYS537 TINYINT,
                  DYS641 TINYINT,
                  DYS472 TINYINT,
                  DYF406S1 TINYINT,
                  DYS511 TINYINT,
                  DYS425 TINYINT,
                  DYS413a TINYINT,
                  DYS413b TINYINT,
                  DYS557 TINYINT,
                  DYS594 TINYINT,
                  DYS436 TINYINT,
                  DYS490 TINYINT,
                  DYS534 TINYINT,
                  DYS450 TINYINT,
                  DYS444 TINYINT,
                  DYS481 TINYINT,
                  DYS520 TINYINT,
                  DYS446 TINYINT,
                  DYS617 TINYINT,
                  DYS568 TINYINT,
                  DYS487 TINYINT,
                  DYS572 TINYINT,
                  DYS640 TINYINT,
                  DYS492 TINYINT,
                  DYS565 TINYINT,
                  DYS710 TINYINT,
                  DYS485 TINYINT,
                  DYS632 TINYINT,
                  DYS495 TINYINT,
                  DYS540 TINYINT,
                  DYS714 TINYINT,
                  DYS716 TINYINT,
                  DYS717 TINYINT,
                  DYS505 TINYINT,
                  DYS556 TINYINT,
                  DYS549 TINYINT,
                  DYS589 TINYINT,
                  DYS522 TINYINT,
                  DYS494 TINYINT,
                  DYS533 TINYINT,
                  DYS636 TINYINT,
                  DYS575 TINYINT,
                  DYS638 TINYINT,
                  DYS462 TINYINT,
                  DYS452 TINYINT,
                  DYS445 TINYINT,
                  YA10 TINYINT,
                  DYS463 TINYINT,
                  DYS441 TINYINT,
                  Y1B07 TINYINT,
                  DYS525 TINYINT,
                  DYS712 TINYINT,
                  DYS593 TINYINT,
                  DYS650 TINYINT,
                  DYS532 TINYINT,
                  DYS715 TINYINT,
                  DYS504 TINYINT,
                  DYS513 TINYINT,
                  DYS561 TINYINT,
                  DYS552 TINYINT,
                  DYS726 TINYINT,
                  DYS635 TINYINT,
                  DYS587 TINYINT,
                  DYS643 TINYINT,
                  DYS497 TINYINT,
                  DYS510 TINYINT,
                  DYS434 TINYINT,
                  DYS461 TINYINT,
                  DYS435 TINYINT)''')

    # calls

    dc.execute('''CREATE TABLE IF NOT EXISTS calls
                  (id LONGINT PRIMARY KEY,
                  person INTEGER REFERENCES people(person),
                  variant INTEGER REFERENCES variants(id),
                  assigned BOOLEAN,
                  positive BOOLEAN,
                  callable BOOLEAN,
                  bq TINYINT,
                  mq TINYINT,
                  nder TINYINT,
                  nanc TINYINT
                  )''')

    # tree

    dc.execute('''CREATE TABLE IF NOT EXISTS tree
                  (id INTEGER PRIMARY KEY,
                  parendid INTEGER,
                  clade CHARACTER(16),
                  variants BLOB,
                  qualities BLOB,
                  ageflags BLOB,
                  children BLOB,
                  ancestralstrs BLOB,
                  originlat REAL,
                  originlong REAL,
                  coverage INTEGER,
                  agecoverage INTEGER,
                  poznikcoverage INTEGER,
                  combbedcoverage INTEGER,
                  olderthan SMALLINT,
                  olderthanunc SMALLINT,
                  olderthankit INTEGER,
                  youngerthan SMALLINT,
                  yongerthanunc SMALLINT,
                  youngerthankit1 INTEGER,
                  youngerthankit2 INTEGER,
                  snpage SMALLINT,
                  snpagelo SMALLINT,
                  snpagehi SMALLINT,
                  snpagepdf BLOB,
                  snpparentpdf BLOB,
                  strage SMALLINT,
                  stragelo SMALLINT,
                  stragehi SMALLINT,
                  stragepdf BLOB,
                  strparentpdf BLOB,
                  combage SMALLINT,
                  combagelo SMALLINT,
                  combagehi SMALLINT,
                  combagepdf BLOB)''')


    # trace (1, "Processing Build 38 BigY files...")
    
def db_inserts(dc,trace,unpack,readVcf):

    # skip to <= 1 - unpack zips

    if (config['skip_to'] <= 1):
        trace (2, "Unpacking ZIP files...")
        #unpack(REDUX_ENV+'/'+config['zip_dir'],REDUX_ENV+'/'+config['unzip_dir'],config['verbosity'])
        unpack()
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
        trace (5, "Associating unzipped files with kits...")
        
    # skip to <= 10 - associate kits with people

    if (config['skip_to'] <= 10):
        trace (2, "Associating kits with people...")
        
    # skip to <= 11 - generate dictionary of variant positions 

    #   The VCF files are parsed twice. The first pass identifies the list
    #   of variants to be queried. The second pass reads the calls for those
    #   variants. This means we only need to treat positions with a variant,
    #   rather than every position in a chromosome.
    #   Using a dictionary here means only one copy of each variant is
    #   created.

    if (config['skip_to'] <= 11):

        #logic 

        trace (2, "Generating database of all variants...")
        vcffiles = [f for f in os.listdir(REDUX_ENV+'/'+config['unzip_dir']) if f.endswith('.vcf')]
        trace (10, "   %i files detected" % len(vcffiles))
        
        variant_dict = {}
        for file in vcffiles:
            vcf_calls = readVcf(REDUX_ENV+'/'+config['unzip_dir']+'/'+ file)
            variant_dict.update(vcf_calls)
        
        trace (10, "   %i variants found" % len(variant_dict))
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)

        # dump variant dictionary into sorted array

        trace (20, "   Dumping variants into array...")
        variant_array = np.array(list(variant_dict.values()))

        # variant_array = np.array([],dtype={'names': ('start', 'anc', 'der'),'formats': ('i4', 'S20', 'S20')})

        trace (30, "      Check variant [0] is %s" % variant_array[0])
        trace (30, "      Check variant [0] position is %s" % variant_array[0][1])
        trace (30, "      Check variant [%s] is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1]))
        trace (30, "      Check variant [%s] position is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1][1]))

        #db calls

        trace (20, "   Inserting data into variant array database...")
        dc.executemany('''INSERT INTO variants(id,ref,alt) VALUES (?,?,?)''', variant_array)

        # Test data has entered database correctly
        # dc.execute('SELECT * FROM variants LIMIT 5')
        # print (dc.fetchone())

        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
        
    # skip to <= 12 - reading calls for variants

    #db calls

    if (config['skip_to'] <= 12):
        trace (2, "Generating database of calls...")
        vcffiles = [f for f in os.listdir(REDUX_ENV+'/'+config['unzip_dir']) if f.endswith('.vcf')]
        trace (10, "   %i files detected" % len(vcffiles))

        dc.execute('''INSERT INTO calls(variant,person)
                      SELECT id, person
                      FROM variants CROSS JOIN people''')

    # Test data has entered database correctly
    # dc.execute('SELECT * FROM calls LIMIT 5')
    # print (dc.fetchone())

    # skip to <= 13 - name variants and derive ancestral values

    # Some variants are positive in the reference sequence, so we need to
    # look up their ancestral values. We'll get the SNP names while we're
    # at it.

    #db calls

    if (config['skip_to'] <= 13):
        trace (2, "Getting names of variants...")

        # Read in SNPs from reference lists
        trace (10, "   Importing SNP reference lists...")

        snp_reference = csv.reader(open(config['b37_snp_file']))
        #for rec in snp_reference:
        #   print "INSERT INTO hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)", (rec[3], rec[4], rec[8], rec[10], rec[11])
        #   dc.execute("INSERT INTO hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)", (rec[3], rec[4], rec[8], rec[10], rec[11]))
        dc.executemany("INSERT INTO hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)",
                       ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
    
        snp_reference = csv.reader(open(config['b38_snp_file']))
        dc.executemany("INSERT INTO hg38(grch38,grch38end,name,anc,der) VALUES (?,?,?,?,?)",
                       ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))

        # Test data has entered database correctly
        #dc.execute('SELECT * FROM hg38 LIMIT 5')
        #print (dc.fetchone())

        # Read in SNPs from reference lists
        # Probably doesn't need done at this point
        # trace (10, "   Joining reference lists to variant database...")


        # dc.execute('''SELECT hg38.grch38, hg38.name
        # FROM hg38
        # INNER JOIN hg19 on hg19.name = hg38.name''')

        # dc.execute('''SELECT variants.id, hg38.name
        # FROM variants
        # LEFT OUTER JOIN hg38 on hg38.grch38 = variants.id''')

        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
            
    # Print final message and exit {{{{

    t = float((time.clock() - start_time))
    trace (1, "Execution finished in: %.3f seconds" % t)

    # }}}

    # sys.exit(0)

