#!/usr/bin/env python3

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

# }}}

REDUX_CONF = os.environ['REDUX_CONF']
REDUX_ENV = os.environ['REDUX_ENV']
config = yaml.load(open(REDUX_CONF))
start_time = 0 # need to fix this

class DB(object):
    
    def __init__(self):
        db = None
        dc = None
        
    def cursor(self,trace):
        trace (1, "Initialising database...")
        self.db = sqlite3.connect('variant.db')
        self.dc = self.db.cursor()
        
    def drop_tables(self):
        if (config['drop_tables'] > 0):
            tbls = "variants hg19 hg38 kits people strs calls tree"
            for tbl in tbls.split():
                self.dc.execute("DROP TABLE IF EXISTS "+tbl)
        
    def create_tables(self):

        #self.dc.execute('create table foo (id int)');
        #self.dc.execute('insert into foo values (1)');
        #sys.exit()

        # variants tbl

        self.dc.execute('''CREATE TABLE IF NOT EXISTS variants
                  (id INTEGER PRIMARY KEY,
                  grch37 INTEGER,
                  length SMALLINT,
                  ref TEXT,
                  alt TEXT,
                  inage BOOLEAN)''')

        # hg19

        self.dc.execute('''CREATE TABLE IF NOT EXISTS hg19
                  (id INTEGER PRIMARY KEY,
                  snp BOOLEAN,
                  grch37 INTEGER,
                  grch37end INTEGER,
                  name CHARACTER(32),
                  anc CHARACTER(64),
                  der CHARACTER(64))''')

        # hg38

        self.dc.execute('''CREATE TABLE IF NOT EXISTS hg38
                  (id INTEGER PRIMARY KEY,
                  snp BOOLEAN,
                  grch38 INTEGER,
                  grch38end INTEGER,
                  name CHARACTER(32),
                  anc CHARACTER(64),
                  der CHARACTER(64))''')

        # kits

        self.dc.execute('''CREATE TABLE IF NOT EXISTS kits
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

        self.dc.execute('''CREATE TABLE IF NOT EXISTS people
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

        self.dc.execute('''CREATE TABLE IF NOT EXISTS strs
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

        self.dc.execute('''CREATE TABLE IF NOT EXISTS calls
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

        self.dc.execute('''CREATE TABLE IF NOT EXISTS tree
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
    
    def insert_variants(self,variant_array):
        self.dc.executemany('''INSERT INTO variants(id,ref,alt) VALUES (?,?,?)''', variant_array)
        # db-debug variants insertion
        #dc.execute('SELECT * FROM variants LIMIT 5')
        #print (dc.fetchone())
        
    def insert_calls(self):
        self.dc.execute('''INSERT INTO calls(variant,person)
            SELECT id, person
            FROM variants CROSS JOIN people''')
        # db-debug calls insertion
        # dc.execute('SELECT * FROM calls LIMIT 5')
        # print (dc.fetchone())
        
    def insert_hg19(self,snp_reference):
        self.dc.executemany("INSERT INTO hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)",
            ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
            
    def insert_hg38(self,snp_reference):
        self.dc.executemany("INSERT INTO hg38(grch38,grch38end,name,anc,der) VALUES (?,?,?,?,?)",
            ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
        # db-debug hg38 insertion
        #self.dc.execute('SELECT * FROM hg38 LIMIT 5')
        #print (dc.fetchone())

