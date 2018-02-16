#!/usr/bin/env python3
#
# Copyright (c) 2018 the Authors
#
# Purpose: Fetch kit data from the Web
#

import yaml
from lib import *

def global_init():
  # environment variable required
  try:
    sys.path.append(os.environ['REDUX_PATH'])
    REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
  except:
    trace(0,"Missing environment variable REDUX_PATH. Aborting.")
    sys.exit()

  # parse the remainder of the configuration settings
  config = yaml.load(open(REDUX_CONF))


def main():
  global_init()
  trace(0, 'Fetching kit info')
  # Get the kits info into json.out no matter what the config says.
  get_kits()
  trace(0, 'Kit info fetched')
  return 0

if __name__ == '__main__':
  sys.exit(main())
