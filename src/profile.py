#!/usr/bin/env python3
# coding: utf-8
#
# Copyright (c) 2018 The Authors
#
# Purpose: profile memory use and timings
#
import psutil
import time, os, yaml
from lib import Trace

REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
config = yaml.load(open(REDUX_CONF))

trace = Trace(level=config['verbosity'])

# Procedure: get_meminfo
# Purpose: return memory usage of this process (resident set size)
def get_meminfo():
    # get process info from OS
    p = psutil.Process(os.getpid())
    # return resident set size
    return p.memory_full_info().rss

# Procedure: profile
# Purpose: decorator for function that adds memory usage and time stat
def profile(myfunc):
    from functools import wraps
    @wraps(myfunc)
    def wrapper(*arg, **kw):
        mem0 = get_meminfo()
        start = time.time()
        result = myfunc(*arg, **kw)
        seconds = time.time() - start
        mem1 = get_meminfo()
        trace(2, '{}: mem: {:,} bytes; time: {:.4f} sec'.format(myfunc.__name__,
                     mem1-mem0, seconds))
        #print('{}: mem: {:,} bytes; time: {:.4f} sec'.format(myfunc.__name__,
        #             mem1-mem0, seconds))
        return result
    return wrapper
