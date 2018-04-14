#!/usr/bin/env python3
# coding: utf-8
#
# Copyright (c) 2018 the Authors
#
# Purpose: routines to maintain and I/O the clade trees
#
import yaml, os
from lib import Trace

REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
config = yaml.load(open(REDUX_CONF))

trace = Trace(level=config['verbosity'])

# Procedure: tree_newclade
# Purpose: create a new clade and fill it with the associated kits and vars
# Return: ID of newly created node
def tree_newclade(dbo, kits, variants, cladename=None):
    if not cladename:
        cladename = str(list(variants)[0])

    dc = dbo.cursor()
    dc.execute('INSERT INTO treeclade(cladeName) VALUES(?)', (cladename,))
    nodeID = dc.lastrowid

    # new node always points to itself at a depth of 0
    dc.execute('''INSERT INTO treepaths (ancestor, descendant, treedepth)
                      values(?,?,0)''', (nodeID, nodeID))
    # INSERT INTO cladekits ...
    tups = []
    for kit in kits:
        tups.append((nodeID, kit))
    dc.executemany('INSERT INTO cladekits(cladeID,pID) values(?,?)', tups)
    # INSERT INTO cladevariants ...
    tups = []
    for var in variants:
        tups.append((nodeID, var))
    dc.executemany('INSERT INTO cladevariants(cladeID,vID) values(?,?)', tups)
    return nodeID

# Procedure: tree_newchild
# Purpose: insert a new child of a specified node
def tree_newchild(dbo, parent, newchild):
    # parent node always points to new child node
    dc = dbo.cursor()
    dc.execute('''INSERT INTO treepaths (ancestor, descendant, treedepth)
                  values(?,?,1)''', (parent, newchild))
    # all parents point to child node
    dc.execute('''INSERT OR IGNORE INTO treepaths
                     (ancestor, descendant, treedepth)
                  SELECT ancestor, ?, treedepth+1 from treepaths
                  where descendant=?''', (newchild, parent))

# Procedure: tree_clade_ancestors
# Purpose: return cladeIDs for all ancestors of a clade
# Return: list of ancestor IDs
def tree_clade_ancestors(dbo, child):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.ancestor
                  WHERE t.descendant=? and t.treedepth > 0''', (child,))
    ids = [c[0] for c in dc]
    return ids

# Procedure: tree_clade_descendants
# Purpose: return cladeIDs for all descendants of a clade
# Return: list of descendant IDs
def tree_clade_descendants(dbo, parent):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.descendant
                  WHERE t.ancestor=? and t.treedepth > 0''', (parent,))
    ids = [c[0] for c in dc]
    return ids

# Procedure: tree_clade_parent
# Purpose: return cladeID for the parent of a clade
# Return: parent ID
def tree_clade_parent(dbo, child):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.ancestor
                  WHERE t.descendant=? and t.treedepth=1''', (child,))
    return dc.fetchone()[0]

# incomplete
# Procedure: tree_matching_clade
def tree_matching_clade(dbo, kitlist, nodelist):
    return
