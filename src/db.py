#!/usr/bin/env python3
# coding: utf-8

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand, Zak Jones
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

import sys,os,sqlite3,yaml,time,csv,json,numpy as np

REDUX_CONF = 'config.yaml'
config = yaml.load(open(REDUX_CONF))


class DB(object):

    def __init__(self, dbfname=config['DB_FILE'], drop=True):
        self.dbfname = dbfname
        # just remove the file, which is often faster than dropping big tables
        if drop and os.path.exists(dbfname):
            os.unlink(dbfname)
        self.db = sqlite3.connect(dbfname)
        self.dc = self.cursor()
        # affect whether or not to wait for data write to disk
        self.dc.execute('PRAGMA synchronous=OFF')
        # force single-user for slightly better performance
        self.dc.execute('PRAGMA locking_mode=EXCLUSIVE')

    def cursor(self):
        return self.db.cursor()

    def run_sql_file(self,FILE):
        with open(FILE,'r') as fh:
            self.dc.executescript(fh.read())

    def commit(self):
        self.db.commit()

    def close(self):
        self.dc.close()

    def create_schema(self, schemafile='schema.sql'):
        self.run_sql_file(os.path.join(config['REDUX_SQL'],schemafile))

    # insert an array of variants into variant definitions table
    # This procedure takes a list or iterator in variant_array, which are
    # pos,anc,der tuples and inserts these into the variants table.  if there
    # are duplicates in variant_array, they are not inserted due to the unique
    # constraint on the variants table.
    def insert_variants(self, variant_array, buildname='hg38'):
        bid = get_build_byname(self.db, buildname)
        alleles = set([v[1] for v in variant_array] + [v[2] for v in variant_array])
        for allele in [a.strip() for a in alleles]:
            db.dc.execute('insert or ignore into alleles(allele) values(?)', (allele,))
        # fixme? may need optimization
        self.dc.executemany('''
            INSERT OR IGNORE INTO variants(buildID,pos,anc,der)
            select b.ID, v.pos, aa.id, ab.id
            from variants v
            cross join build b on b.buildNm=?
            cross join allele ab on ab.allele=?
            cross join allele aa on aa.allele=?''',
            [(bid,)+v for v in variant_array])

    # insert a vector of variant ids to insert for a given person specified by
    # pID. This procedure inserts ids into the "calls" table and operates with
    # variant ids (i.e. you already have the variant ids)
    def insert_calls(self, pid, calls):
        self.dc.executemany('INSERT INTO vcfcalls(pID,vID) values (?,?)',
                                [(pid,v) for v in calls])




# test framework
if __name__=='__main__':
    db = DB()
    db.create_schema()
    db.commit()
