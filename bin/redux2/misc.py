#!/bin/python

# authors/licensing {{{

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import os

# }}}

REDUX_CONF = os.environ['REDUX_CONF']
REDUX_ENV = os.environ['REDUX_ENV']

# routines - debug

def trace (level, msg):
    print(msg)
    #if level <= config['verbosity']:
    #    print(msg)
    
