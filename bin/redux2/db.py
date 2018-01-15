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
REDUX_SQL = os.environ['REDUX_SQL']
config = yaml.load(open(REDUX_CONF))
start_time = 0 # need to fix this

class DB(object):
    
    def __init__(self):
        db = None
        dc = None
        
    def cursor(self):
        #trace (1, "Initialising database...")
        self.db = sqlite3.connect('variant.db')
        self.dc = self.db.cursor()
        
    def run_sql_file(self,FILE):
        fh = open(REDUX_SQL+'/'+FILE,'r');
        try:
            #print(fh.read())
            self.dc.executescript(fh.read())
        finally:
            fh.close()
    def commit(self):
        self.db.commit()

    #redux2 ddl+dml

    def drop_v1_tables(self):
        self.run_sql_file('redux2-drop-tables.sql')
        
    def create_v1_tables(self):
        self.run_sql_file('redux2-schema.sql')
        
    def insert_v1_variants(self,variant_array):
        self.dc.executemany('''INSERT INTO v1_variants(id,ref,alt) VALUES (?,?,?)''', variant_array)
        # db-debug variants insertion
        #dc.execute('SELECT * FROM v1_variants LIMIT 5')
        #print (dc.fetchone())
        
    def insert_v1_calls(self):
        self.dc.execute('''INSERT INTO v1_calls(variant,person)
            SELECT id, person
            FROM v1_variants CROSS JOIN v1_people''')
        # db-debug calls insertion
        # dc.execute('SELECT * FROM v1_calls LIMIT 5')
        # print (dc.fetchone())
        
    def insert_v1_hg19(self,snp_reference):
        self.dc.executemany("INSERT INTO v1_hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)",
            ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
            
    def insert_v1_hg38(self,snp_reference):
        self.dc.executemany("INSERT INTO v1_hg38(grch38,grch38end,name,anc,der) VALUES (?,?,?,?,?)",
            ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
        # db-debug hg38 insertion
        #self.dc.execute('SELECT * FROM v1_hg38 LIMIT 5')
        #print (dc.fetchone())

    #v2 db schema ddl+dml

    def drop_tables(self):
        self.run_sql_file('v2-drop-tables.sql')
        
    def create_tables(self):
        self.run_sql_file('v2-schema.sql')

    #tree sort prototype ddl+dml
    
    def drop_sort_tables(self):
        self.run_sql_file('sort-drop-tables.sql')
        
    def create_sort_tables(self):
        self.run_sql_file('sort-schema.sql')
        
    def insert_sort_data(self):
        foo = 2

