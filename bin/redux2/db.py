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

import sys,os,sqlite3,yaml,time,csv,json,numpy as np

# }}}

REDUX_CONF = os.environ['REDUX_CONF']
config = yaml.load(open(REDUX_CONF))
start_time = 0 # need to fix this

#TODO: need to put this info properly into yaml file (for now, I'm hacking bashrc)
REDUX_ENV = os.environ['REDUX_ENV']
REDUX_SQL = os.environ['REDUX_SQL']
REDUX_DATA = os.environ['REDUX_DATA']

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
        
    def insert_sample_sort_data(self):

        #sample data: 3019783,M343,1,1,Null,1,1,1,1,1,1,1
        #{'v': '3019783', 'n': 'M343', 'k1': '1', 'k2': '1', 'k3': 'Null', 'k4': '1', 'k5': '1', 'k6': '1', 'k7': '1', 'k8': '1', 'k9': '1', 'k10': '1', 'k11': None, 'k12': None}

        #create table s_kits(
        # kit_id  int,  -- later this can be person_id
        # sort_order int
        #);

        #create table s_variants (
        # -- variant_id int, -- not needed for prototype
        # variant_loc int,  -- PK
        # name varchar(20)
        # -- old_reference varchar(2), -- commenting out right now cuz not part of ian's doc
        #);

        #create table s_calls(
        # kit_id int,
        # variant_loc int,
        # assigned boolean
        #);

        cols=10
        for k in range(1,cols+1):
            self.dc.execute("insert into s_kits (kit_id) values ("+str(k)+");")

        with open(REDUX_DATA+'/sample-sort-data.csv','r') as FILE:
            for row in csv.DictReader(FILE,'v n k1 k2 k3 k4 k5 k6 k7 k8 k9 k10'.split()):
                row = json.loads(json.dumps(row).replace('\\ufeff','')) #hack: remove byte order mark
                self.dc.execute("insert into s_variants (variant_loc,name) values ("+row['v']+",'"+row['n']+"');")
                #print(' - inserting sample variant data: '+str(row['v']))
                for k in range(1,cols+1):
                    kv = str(row['k'+str(k)])
                    #'null' if kv == "None" else kv
                    vv = str(row['v'])
                    #print (kv)
                    self.dc.execute("insert into s_calls (kit_id,variant_loc,assigned) values ("+str(k)+","+vv+","+kv+");")
                    #self.commit()
                    #print (kv+":"+vv)
                #break;
                #sys.exit()
                #print(row)
                #print(row.encode('utf-8-sig'))
            #for (l,n,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12) in dr]
            #    print (n)
            #self.dc.executemany("insert into s_kits (col1, col2) VALUES (?, ?);", to_db)
            #(variant_loc,name,) = t_db
            #con.commit()
            #con.close()
        self.commit()
        
    def sort_data(self):

        #sql_1 = "select kit_id,count(*) as pos_k_cnt from s_calls where assigned = 1 group by kit_id order by count(*) desc;"
        #self.dc.execute(sql_1)
        #kitA = self.dc.fetchall()
        #print("---")
        #print("kit check")
        #print(kitA)
        #[(4, 9), (9, 8), (2, 7), (7, 7), (8, 7), (3, 6), ... ]

        print("===")
        print("START: 'For Each Variant'")
        print("===")
        sql_0 = "select distinct variant_loc from s_calls;"
        self.dc.execute(sql_0)
        A = self.dc.fetchall()
        print("A - distinct variants")
        print(A)
        print("---")
        print('numA - num distinct variants')
        print(len(A)) #num distinct variants

        print("===")
        print("SORT STEP: B")
        print("===")
        sql_1 = "select distinct variant_loc from s_calls where assigned = 0;"
        self.dc.execute(sql_1)
        B1 = self.dc.fetchall()
        B0 = set(A)-set(B1)
        #print("---")
        #print("sql result")
        #print(f1)
        #sys.exit()
        #print("---")
        print("B0 - variants that don't have negs")
        #print(B0) #B0 -  variants that don't have negs 
        print(list(B0))
        print("---")
        print("B1 - variants that have negs")
        print(B1) #B1 - variants that have negs 

        print("===")
        print("SORT STEP: C")
        print("===")
        sql_2 = "select variant_loc,count(*) as cnt from s_calls where assigned = 1 group by variant_loc;"
        self.dc.execute(sql_2)
        f2 = self.dc.fetchall()
        f2a = list(filter(lambda x: x[1]==(len(f2)-1), f2))
        C1 = list(set(B1) & set(f2a)) # intersection = singletons
        C0 = list(set(B1)-set(C1))
        print("list of *all* one person +ve's")
        print(f2a)
        print("---")
        print("not singletons")
        print(C0)
        print("---")
        print("singletons")
        print(C1)

        print("===")
        print("SORT STEP: D")
        print("===")
        sql_3 = "select distinct variant_loc from s_calls where assigned is null group by variant_loc;"
        self.dc.execute(sql_3)
        f3 = self.dc.fetchall()
        D0 = list(set(C1)-set(f3)) #imperfect variants
        D1 = list(set(f3)-set(D0))
        print("list of variants that are sometimes not called")
        print(f3)
        print("---")
        print("imperfect variants")
        print(D0)
        print("---")
        print("calls of perfect share variants")
        print(D1)

        sys.exit()

        #sql_2b = "select variant_loc,count(*) as pos_v_cnt from s_calls where assigned = 0 group by variant_loc order by count(*) desc;"
        #self.dc.execute(sql_2b)
        #varAn = self.dc.fetchall()
        #print("---")
        #print("variant negative check")
        #print(varAn)

        sql_2b = "select variant_loc,count(*) as pos_v_cnt from s_calls where assigned = 0 group by variant_loc order by count(*) desc;"
        sql_2c = "select variant_loc,count(*) as pos_v_cnt from s_calls where assigned is not null group by variant_loc order by count(*) desc;"
        self.dc.execute(sql_2c)
        varAa = self.dc.fetchall()
        print("---")
        print("variant all check")
        print(varAa)
        #(3) 9 perfectly called variants - execute sort on these
        #(4) 6 imperfectly called variants - do Step A

        sql_3 = "select * from s_calls order by kit_id,assigned;"
        self.dc.execute(sql_3)
        callsA = self.dc.fetchall()
        print("---")
        #[(1, 12060401, None), (1, 6920349, 0), (1, 7378685, 0), (1, 13668461, 0), (1, 19538924, 0), ... ]
        print (callsA);

        #Note: build the default structure with the kits ordered like kitA and the variants ordered like varA
        #[{"k1":[{"v12060401",1)},{"v6920349",1), ... ]}
        #[{"k2":[{"v12060401",1)},{"v6920349",None), ... ]}
        #and display it
        #for call in calls:
        #   for K in kits:
        #    ...
        #   sort_positive_variants(kit_id)



