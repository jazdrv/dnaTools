#!/usr/bin/env python3

# license + libs{{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

import sys,os,sqlite3,yaml,csv,json,numpy as np

# }}}

try:
    config = yaml.load(open(os.environ['REDUX_CONF_ZAK']))
except:
    print("Missing environment variable REDUX_CONF_ZAK. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

class DB1(object):
    
    def __init__(self):
        self.db = None
        self.dc = None
        
    def db_init(self):
        return sqlite3.connect(config['DB_FILE'])
        
    def cursor(self):
        return self.db.cursor()
        
    def sql_exec_file(self,FILE):
        fh = open(config['REDUX_SQL']+'/'+FILE,'r');
        try:
            self.dc.executescript(fh.read())
        finally:
            fh.close()
        
    def sql_exec(self,sql,dc=None):
        if config['DBG_SQL']: print("[SQL] "+sql)
        if dc is not None:
            dc.execute(sql)
        else:
            self.dc.execute(sql)
            if config['COMMITS_ON']: self.commit()
        
    def sql_exec_many(self,sql,itms,dc=None):
        if config['DBG_SQL']: print("[SQL] "+sql)
        if dc is not None:
            dc.executemany(sql,itms)
        else:
            self.dc.executemany(sql,itms)
            if config['COMMITS_ON']: self.commit()
        
    def fetchall(self):
        return self.dc.fetchall()
        
    def fetchone(self):
        return self.dc.fetchone()
        
    def commit(self):
        self.db.commit()
