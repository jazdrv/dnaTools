import unittest,sys
from lib import *
from db import *

class TestDB(unittest.TestCase):

    def test_db(self):
        cur = db_init(trace)
        db_drop_tables(cur)
        db_create_tables(cur)
        db_inserts(cur,trace,unpack,readVcf)
        #res = fact(5)
        #self.assertEqual(res, 120)

if __name__ == '__main__':
    unittest.main()
