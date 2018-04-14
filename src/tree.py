#!/usr/bin/env python3
# coding: utf-8
#
# Copyright (c) 2018 the Authors
#
# Purpose: routines to maintain and I/O the clade trees
#
import yaml, os
from lib import Trace
from array_api import *

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

# Procedure: tree_node_kits
# Purpose: return list of kits associated with node
def tree_node_kits(dbo, node):
    dc = dbo.cursor()
    dc.execute('select pid from cladekits where cladeid=?', (node,))
    pidlist = []
    for (pid,) in dc:
        pidlist.append(pid)
    return get_kit_ids(dbo, pidlist).values()

# Procedure: tree_node_variants
# Purpose: return list of kits associated with node
def tree_node_variants(dbo, node):
    dc = dbo.cursor()
    dc.execute('select vid from cladevariants where cladeid=?', (node,))
    vidlist = []
    for (vid,) in dc:
        vidlist.append(vid)
    return get_variant_snpnames(dbo, vidlist).values()


# Procedure: tree_to_dot
# Purpose: turn a tree into DOT language
def tree_to_dot(dbo, top):
    def nodestr(nodeid, kits, variants):
        return '{} [label="{}"];\n'.format(str(nodeid), '\\n'.join([str(v) for v in variants if v]))
    def linkstr(parentid, childid):
        return '{} -> {};\n'.format(str(parentid), str(childid))
    nodelist = tree_clade_descendants(dbo, top)
    out = ''
    for node in nodelist:
        out += nodestr(node, tree_node_kits(dbo,node),
                           tree_node_variants(dbo, node))
    dc = dbo.cursor()
    dc.execute('select ancestor, descendant from treepaths where treedepth=1')
    pairs = [(r[0],r[1]) for r in dc if r[0] in nodelist or r[1] in nodelist]
    for pair in pairs:
        out += linkstr(pair[0], pair[1])
    return '''digraph tree {
edge [style=bold,color=blue];
node [fontname="Helvetica" fillcolor=white shape=none margin="0,0"];

''' + out + '\n}'
