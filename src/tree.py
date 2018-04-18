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
from profile import profile

REDUX_CONF = os.path.join(os.environ['REDUX_PATH'], 'config.yaml')
config = yaml.load(open(REDUX_CONF))

trace = Trace(level=config['verbosity'])

# Procedure: tree_get_treetops
# Purpose: return all tree nodes labeled "Top"
def tree_get_treetops(dbo):
    dc = dbo.cursor()
    dc.execute('select id from treeclade where cladename=?', ('Top',))
    return [ID[0] for ID in dc]

# Procedure: tree_delete_tree
# Purpose: remove all remnants of a tree starting at treenode
def tree_delete_tree(dbo, treenode):
    dc = dbo.cursor()
    dc.execute('''delete from treeclade where id in
                  (select descendant from treepaths where ancestor=?)''',
                  (treenode,))
    dc.execute('''delete from cladevariants where cladeid in
                  (select descendant from treepaths where ancestor=?)''',
                  (treenode,))
    dc.execute('''delete from cladekits where cladeid in
                  (select descendant from treepaths where ancestor=?)''',
                  (treenode,))
    dc.execute('''delete from treepaths where treedepth=0 and ancestor not in
                  (select id from treeclade)''')

    return

# Procedure: tree_newclade
# Purpose: create a new clade and fill it with the associated kits and vars
# Return: ID of newly created node
def tree_newclade(dbo, kits, variants, cladename=None):
    if not cladename:
        try:
            cladename = str(list(variants)[0])
        except IndexError:
            trace(0, 'Unknown variant: {},{}'.format(kits,variants))
            cladename='Unknown'

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
    dbo.commit()
    return nodeID

# Procedure: tree_clade_update_variants
# Purpose: replace clade's variants with a different set of variants
def tree_clade_update_variants(dbo, clade, variants):
    dc = dbo.cursor()
    dc.execute('delete from cladevariants where cladeid=?', clade)
    newlist = [(cladeid, v) for v in variants]
    dc.executemany('insert into cladevariants values(?,?)', newlist)
    return

# Procedure: tree_add_child_clade
# Purpose: insert a new child of a specified node
def tree_add_child_clade(dbo, parent, newchild):
    # parent node always points to new child node
    dc = dbo.cursor()
    dc.execute('''INSERT INTO treepaths (ancestor, descendant, treedepth)
                  values(?,?,1)''', (parent, newchild))
    # all parents point to child node
    dc.execute('''INSERT OR IGNORE INTO treepaths
                     (ancestor, descendant, treedepth)
                  SELECT ancestor, ?, treedepth+1 from treepaths
                  where descendant=?''', (newchild, parent))
    dbo.commit()
    # child clade contains variants of parent - remove them from parent
    pvars = set(tree_clade_get_variants(dbo, parent))
    if not pvars:
        trace(0, 'NO VARIANTS: parent node {}'.format(parent))
        return
    cvars = set(tree_clade_get_variants(dbo, newchild))
    newvars = pvars - cvars
    if not newvars or newvars == pvars:
        trace(0, 'ERROR: child {} not in {}'.format(newchild, parent))
        trace(1, 'cvars: {}'.format(cvars))
        trace(1, 'pvars: {}'.format(pvars))
        trace(1, 'newvars: {}'.format(newvars))

    tree_clade_update_variants(dbo, parent, newvars)

# Procedure: tree_clade_delete_node
# Purpose: delete a leaf node
def tree_clade_delete_node(dbo, oldnode):
    # parent node always points to new child node
    dc = dbo.cursor()
    dc.execute('''INSERT INTO treepaths (ancestor, descendant, treedepth)
                  values(?,?,1)''', (parent, newchild))
    # all parents point to child node
    dc.execute('''INSERT OR IGNORE INTO treepaths
                     (ancestor, descendant, treedepth)
                  SELECT ancestor, ?, treedepth+1 from treepaths
                  where descendant=?''', (newchild, parent))
    dbo.commit()

# Procedure: tree_clade_get_ancestors
# Purpose: return cladeIDs for all ancestors of a clade
# Return: list of ancestor IDs
def tree_clade_get_ancestors(dbo, child):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.ancestor
                  WHERE t.descendant=? and t.treedepth > 0''', (child,))
    ids = [c[0] for c in dc]
    return ids

# Procedure: tree_clade_get_descendants
# Purpose: return cladeIDs for all descendants of a clade
# Return: list of descendant IDs
def tree_clade_get_descendants(dbo, parent):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.descendant
                  WHERE t.ancestor=? and t.treedepth > 0''', (parent,))
    ids = [c[0] for c in dc]
    return ids

# Procedure: tree_clade_get_children
# Purpose: return cladeIDs for all direct children of a clade
# Return: list of child IDs
def tree_clade_get_children(dbo, parent):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.descendant
                  WHERE t.ancestor=? and t.treedepth=1''', (parent,))
    ids = [c[0] for c in dc]
    return ids

# Procedure: tree_clade_get_parent
# Purpose: return cladeID for the parent of a clade
# Return: parent ID
def tree_clade_get_parent(dbo, child):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.ancestor
                  WHERE t.descendant=? and t.treedepth=1''', (child,))
    return dc.fetchone()[0]

# Procedure: tree_insert_clade
# Purpose: insert a new clade between parent and child
@profile
def tree_insert_clade(dbo, newclade, parent, child):
    tree_add_child_clade(dbo, parent, newclade)
    dc = dbo.cursor()
    # remove old paths to child, add child under newclade
    dc.execute('delete from treepaths where descendant=? and treedepth>0',
                   (child,))
    tree_add_child_clade(dbo, newclade, child)
    # all children of previous child are now depth+1
    dc.execute('''update treepaths set treedepth=treedepth+1
                  where descendant in
                  (select descendant from treepaths
                     where ancestor=? and treedepth > 0)''', (child,))
    # all children of previous child are also a descendant of newclade
    dc.execute('''insert into treepaths(ancestor,descendant,treedepth)
                  select ?,descendant,treedepth+1 from treepaths
                  where ancestor=? and treedepth>0''', (newclade, child))
    # add child kits into new clade and add newclade kits into parent
    tree_clade_add_kits(dbo, newclade, tree_clade_get_kits(dbo,child))
    tree_clade_add_kits(dbo, parent, tree_clade_get_kits(dbo,newclade))

# Procedure: tree_clade_get_kits
# Purpose: return list of kits associated with node
def tree_clade_get_kits(dbo, node):
    dc = dbo.cursor()
    dc.execute('select pid from cladekits where cladeid=?', (node,))
    pidlist = []
    for (pid,) in dc:
        pidlist.append(pid)
    return pidlist

# Procedure: tree_clade_get_kit_defs
# Purpose: return dict of kits definitions associated with node
def tree_clade_get_kit_defs(dbo, node):
    dc = dbo.cursor()
    dc.execute('select pid from cladekits where cladeid=?', (node,))
    pidlist = []
    for (pid,) in dc:
        pidlist.append(pid)
    return get_kit_ids(dbo, pidlist)

# Procedure: tree_clade_add_kits
# Purpose: add a list of kits to a given node
def tree_clade_add_kits(dbo, node, newkits):
    nodekits = tree_clade_get_kits(dbo, node)
    allkits = set(nodekits + newkits)
    dc = dbo.cursor()
    # rollback update if sanity check fails
    dbo.commit()
    # also push kits up the tree to all parent nodes
    for ancestor in tree_clade_get_ancestors(dbo,node) + [node]:
        trace(3, 'allkits: {}'.format(allkits))
        for pid in allkits:
            dc.execute('''insert or ignore into cladekits(cladeid,pid)
                          values(?,?)''', (node,pid))
    # FIXME - check here for any inconsistent calls up the tree
    # roll-back update if inconsistencies found
    dbo.commit()
    return

# Procedure: tree_clade_get_variants
# Purpose: return list of variant ids associated with node
def tree_clade_get_variants(dbo, node):
    dc = dbo.cursor()
    dc.execute('select vid from cladevariants where cladeid=?', (node,))
    vidlist = []
    for (vid,) in dc:
        vidlist.append(vid)
    return vidlist

# Procedure: tree_clade_get_variant_defs
# Purpose: return dict of variant definitions associated with node
def tree_clade_get_variant_defs(dbo, node):
    dc = dbo.cursor()
    dc.execute('select vid from cladevariants where cladeid=?', (node,))
    vidlist = []
    for (vid,) in dc:
        vidlist.append(vid)
    trace(3, 'node: {}; vidlist: {}'.format(node, vidlist))
    return get_variant_defs(dbo, vidlist)

# Procedure: tree_is_leaf_node
# Purpose: identify a leaf node (has no children)
def tree_is_leaf_node(dbo, node):
    if tree_clade_get_children(dbo, node):
        return False
    return True

# Procedure: tree_clades_depth_first
# Purpose: return list of tree clade ids, deepest first
def tree_clades_depth_first(dbo, node):
    dc = dbo.cursor()
    dc.execute('''select descendant from treepaths where ancestor=?
                  order by treedepth desc''', (node,))
    return list([n[0] for n in dc])

# Procedure: tree_clade_delete_node
# Purpose: delete a leaf node
def tree_clade_delete_node(dbo, node):
    if not tree_is_leaf_node(node):
        return
    dc = dbo.cursor()
    dc.execute('delete from treepaths where descendant=?', (node,))
    dc.execute('delete from treeclade where id=?', (node,))
    dc.execute('delete from cladevariants where cladeid=?', (node,))
    dc.execute('delete from cladekits where cladeid=?', (node,))
    dbo.commit()
    return

# Procedure: tree_func_bottom_up
# Purpose: apply a function to each node in the tree, bottom up (leaf first)
def tree_func_bottom_up(dbo, node, func):
    dc = dbo.cursor()
    if tree_is_leaf_node(dbo, node):
        func(node)
        return
    else:
        for child in tree_clade_get_children(dbo, node):
            tree_func_bottom_up(dbo, child, func)
    func(node)

# Procedure: tree_merge_into
# Purpose: merge one tree into another tree
# Info:
#   The "into" tree should be the top node of a tree to receive the new tree
#   data. The "from" tree should be a node that we will attempt to merge into
#   the receiving tree. When a node can't be merged, it's added to the
#   remainder tree.
@profile
def tree_merge_into(dbo, fromtree, intotree, remainder):
    # if the "from" tree doesn't have any children, merge it and return
    if tree_is_leaf_node(dbo, fromtree):
        tree_clade_merge_into(dbo, fromtree, intotree, remainder)
        return

    # process all of the children of the "from" tree
    for child in tree_clade_get_children (dbo, fromtree):
        trace(3, 'merge child: {}'.format(child))
        if tree_is_leaf_node(dbo, child):
            tree_clade_merge_into(dbo, child, intotree, remainder)
        else:
            tree_merge_into(dbo, child, intotree, remainder)

    # after processing all of this node's children, process this node
    tree_clade_merge_into(dbo, fromtree, intotree, remainder)
    return

# Procedure: tree_clade_merge_into
# Purpose: merge a node into another tree
# Info:
#   WORK IN PROGRESS
#   This procedure is currently mostly rubbish and needs much more work. A
#   clade is defined as a set of variants that a set of kits has. To merge it
#   into another tree, the variants in common with some node of the other tree
#   need to be factored out, then the node needs to be merged under that node
#   that had some of the variants. If no variants are found in common anywhere
#   in the tree, the new node needs to be parented under the top node of the
#   tree. When some variants are shared with the parent, the parent needs to be
#   split into two, with the new parent being all of the shared variants.
@profile
def tree_clade_merge_into(dbo, fromnode, intotree, remainder):

    trace(2, 'merge {} into {}'.format(fromnode, intotree))

    fromvars = set(tree_clade_get_variants(dbo, fromnode))
    fromkits = tree_clade_get_kits(dbo, fromnode)

    # skip over treetop node or other invalid nodes
    if not fromkits or not fromvars:
        trace(2, 'empty node skipped: {},{}'.format(fromnode,intotree))
        return

    # traverse the target tree from the bottom up (leaf nodes first)
    for node in tree_clades_depth_first(dbo, intotree):
        intovars = set(tree_clade_get_variants(dbo, node))

        # three cases: 1) set of variants in "from" node identically matches
        # the set of variants of some node in "to" tree; 2) "from" set is a
        # subset of some node's variants in "to" tree; 3) there are no shared
        # variants with any node in "to" tree. NB: the algorithm assumes that
        # we already handled all nodes below (descendants of) fromnode, but
        # there's still a complication in that fromnode might represent a span
        # of nodes in the "to" tree. E.g. if "to" tree has a>b>c and "from"
        # node has a,b,c, we should still merge the from node. This might best
        # be handled as case 4) "from" node has a superset of the variants in a
        # node of "to" tree. Then we can look up the tree at the ancestors and
        # ensure the variants are above and merge if yes.

        # CASE 1: variants are identical to some node - copy kits and be done
        if intovars == fromvars:
            tree_clade_add_kits(dbo, node, fromkits)
            return

        # Came to a node that contains all new variants.
        # No child of this node contains all since we went depth first.
        # NOTE: should probably handle case #4 above this one
        elif fromvars.issubset(intovars):
            allchildvars = set()
            for child in tree_clade_get_children(dbo, node):
                variants = set(tree_clade_get_variants(dbo, child))
                if fromvars.issuperset(variants):
                    # all of the child's variants are in new node
                    # in hierarchy, only the difference vars are in newclade
                    newvars = fromvars - variants
                    newclade = tree_newclade(dbo, fromkits, list(newvars))
                    # insert this new clade between the node and its child
                    tree_insert_clade(dbo, newclade, node, child)
                    return
                else:
                    allchildvars = allchildvars.union(variants)
            # was not a superset for any child - add as child of node
            if allchildvars.intersection(fromvars):
                trace(0, 'INCONSISTENT new,parent,vars {},{},{}'.format(
                    fromnode, node, allchildvars))
            newclade = tree_newclade(dbo, fromkits, list(fromvars))
            tree_add_child_clade(dbo, node, newclade)
            tree_clade_add_kits(dbo, node, fromkits)
            return

    # reaching here means no parent node found in intotree
    # fromnode might be a parent node of one of the treetop nodes
    for node in tree_clade_get_children(dbo, intotree):
        intovars = set(tree_clade_get_variants(dbo, node))
        if intovars.issubset(fromvars):
            # FIXME add intotree as parent of fromnode
            # and reparent node as child of fromnode
            # and add kits of node to fromnode
            # and check for inconsistent calls
            pass
    dbo.commit()

# Procedure: tree_to_dot
# Purpose: turn a tree into DOT language
def tree_to_dot(dbo, top):
    def nodestr(nodeid, kits, variants):
        return '{} [label="{}"];\n'.format(str(nodeid), '\\n'.join([str(v) for v in variants]+['kit'+str(k) for k in kits]))
    def linkstr(parentid, childid):
        return '{} -> {};\n'.format(str(parentid), str(childid))
    nodelist = tree_clade_get_descendants(dbo, top)
    out = ''
    for node in nodelist:
        out += nodestr(node, list(tree_clade_get_kit_defs(dbo,node).values()),
               [v[3] for v in tree_clade_get_variant_defs(dbo, node).values()])
    dc = dbo.cursor()
    dc.execute('select ancestor, descendant from treepaths where treedepth=1')
    pairs = [(r[0],r[1]) for r in dc if r[0] in nodelist or r[1] in nodelist]
    for pair in pairs:
        out += linkstr(pair[0], pair[1])
    return '''digraph tree {
edge [style=bold,color=blue];
node [fontname="Helvetica" fillcolor=white shape=none margin="0,0"];

''' + out + '\n}'

# test framework
if __name__=='__main__':
    from db import DB
    dbo = DB(drop=False)
    for tree in tree_get_treetops(dbo):
        tree_delete_tree(dbo, tree)
    dbo.commit()
    sys.exit(0)

    tree_merge_into(dbo,149,1,None)
    tree_merge_into(dbo,302,1,None)
    tree_merge_into(dbo,447,1,None)
    tree_merge_into(dbo,590,1,None)
    tree_merge_into(dbo,740,1,None)
    tree_merge_into(dbo,888,1,None)
    tree_merge_into(dbo,1035,1,None)
    tree_merge_into(dbo,1190,1,None)
    tree_merge_into(dbo,1343,1,None)
    tree_merge_into(dbo,1498,1,None)
    tree_merge_into(dbo,1647,1,None)
    tree_merge_into(dbo,1797,1,None)
    tree_merge_into(dbo,1950,1,None)
    tree_merge_into(dbo,2097,1,None)
    tree_merge_into(dbo,2251,1,None)
    tree_merge_into(dbo,2398,1,None)
    tree_merge_into(dbo,2552,1,None)
    tree_merge_into(dbo,2703,1,None)
    tree_merge_into(dbo,2859,1,None)
    tree_merge_into(dbo,3008,1,None)
    tree_merge_into(dbo,3159,1,None)
    tree_merge_into(dbo,3240,1,None)
    tree_merge_into(dbo,3396,1,None)
    tree_merge_into(dbo,3550,1,None)
    trace(0, 'writing tree to tree.dot...')
    open('tree.dot','w').write(tree_to_dot(dbo, 1))
    sys.exit(0)

    treetop = tree_newclade(dbo, [], [], 'Top')
    c1 = tree_newclade(dbo, [1,2,3,4], [5,6,7,8])
    tree_add_child_clade(dbo, treetop, c1)
    c2 = tree_newclade(dbo, [1,2], [5,6])
    tree_add_child_clade(dbo, c1, c2)

    treetop2 = tree_newclade(dbo, [], [], 'Top')
    c1 = tree_newclade(dbo, [30,31,32], [5,6,7,8])
    tree_add_child_clade(dbo, treetop2, c1)
    c2 = tree_newclade(dbo, [30,31], [5,6])
    tree_add_child_clade(dbo, c1, c2)
    c3 = tree_newclade(dbo, [10,20,40,41], [5,6,7])
    tree_insert_clade(dbo, c3, c1, c2)

    tree_merge_into(dbo,treetop2,treetop,None)

    open('tree.dot','w').write(tree_to_dot(dbo, treetop))


    tree_func_bottom_up(dbo, 1, print)
    print(tree_clades_depth_first(dbo,1))
    print(tree_clade_get_ancestors(dbo,4))
