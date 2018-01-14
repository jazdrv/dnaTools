#!/usr/bin/env python3
'''
Purpose: use the API for for DNA data warehouse to get kit metadata and
populate database

For free distribution under the terms of the GNU General Public License,
version 3 (29 June 2007)
https://www.gnu.org/licenses/gpl.html

Jef Treece, 11 Jan 2018
'''
import json, requests, sqlite3


# pull information about the kits from the web api
def get_kits (API='http://haplogroup-r.org/api/v1/uploads.php', qry='format=json'):
    url = '?'.join([API, qry])
    try:
        res = requests.get(url)
        js = json.loads(res.content)
    except:
        print('Failed to pull kit metadata from {}'.format(API))
        raise # fixme - what to do on error
    return js

# update the information about the kits
def update_metadata(db, js):
    rows = []
    for kit in js:
        rows.append((kit['kitId'], kit['surname'], kit['country']))
    dc = db.cursor()
    # fixme - handle updating existing data
    dc.executemany('insert or ignore into people(kitid,mdkaname,mdkacountry) values (?,?,?)', rows)
    dc.close()

def create_metadata_table(db):
    dc = db.cursor()
    dc.execute('''create table if not exists people
                  (person INTEGER PRIMARY KEY,
                  kitid TEXT,
                  mdkaname TEXT,
                  mdkacountry TEXT)''')
    dc.close()

def dump_db(db):
    dc = db.cursor().execute('select * from people order by mdkaname')
    for id,kid,sn,co in dc:
        if not sn:
            sn="None"
        if not co:
            co="None"
        print('{:5} = {:20}{:20}{:20}'.format(id,kid,sn,co))

def main():
    db = sqlite3.connect('metadata.db')
    create_metadata_table(db)
    js = get_kits()
    print('Number of kits: {}'.format(len(js)))
    update_metadata(db, js)
    db.commit()
    dump_db(db)


main()
