#!/usr/bin/env python3
# coding: utf-8
#
# Copyright (c) 2018 the Authors
#
# Purpose: low level interface to the data layer for DNA analysis
#

import os, sqlite3, yaml

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


# test framework
if __name__=='__main__':
    db = DB()
    db.create_schema()
    db.commit()
