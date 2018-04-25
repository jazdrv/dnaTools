#!/usr/bin/env python3
# coding: utf-8
#
# Copyright (c) 2018 the Authors
#
# Purpose: routines to maintain and I/O clade trees
# Clade trees are a hierarchical tree with kitid,variant ids at each node
# Practical examples in unit tests at end of this file
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
def tree_get_treetops(dbo, node=None):
    dc = dbo.cursor()
    if not node:
        dc.execute('select id from treeclade where cladename=?', ('Top',))
    else:
        dc.execute('''select distinct ancestor from treepaths p
                      inner join treeclade c on c.id=p.ancestor
                      where p.descendant=? and c.cladename=?''',
                      (node, 'Top'))
    return [ID[0] for ID in dc]

# Procedure: tree_cleanup
# Purpose: clean up orphaned entries from the tree tables
def tree_cleanup(dbo):
    dc = dbo.cursor()
    dc.execute('''delete from treepaths where ancestor not in
                  (select id from treeclade)''')
    dc.execute('''delete from cladekits where cladeid not in
                  (select id from treeclade)''')
    dc.execute('''delete from cladevariants where cladeid not in
                  (select id from treeclade)''')
    dc.execute('''delete from treeclade where id not in
                  (select ancestor from treepaths
                       union
                   select descendant from treepaths)''')
    dbo.commit()


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
    dc.execute('''delete from treepaths where ancestor=? or descendant=?''',
                  (treenode,treenode))
    dbo.commit()
    # there will still be some nodes pointing to themselves in treepaths
    tree_cleanup(dbo)
    return

# Procedure: tree_clone_tree
# Purpose: copy the subtree starting at fromtree under the destination node
def tree_clone_tree(dbo, fromtree, destination):
    nvars = tree_clade_get_variants(dbo, fromtree)
    nkits = tree_clade_get_kits(dbo, fromtree)
    newchild = tree_newclade(dbo, nkits, nvars)
    tree_add_child_clade(dbo, destination, newchild)
    for child in tree_clade_get_children(dbo, fromtree):
        tree_clone_tree(dbo, child, newchild)
    dbo.commit()
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
    trace(9, 'tree_newclade created node {} for kits {} and variants {}'.format(
        nodeID, kits, variants))
    dbo.commit()
    return nodeID

# Procedure: tree_node_depth
# Purpose: return max depth of this node within any tree
def tree_node_depth(dbo, node):
    dc = dbo.cursor()
    dc.execute('select max(treedepth) from treepaths where descendant=?',
                   (node,))
    return dc.fetchone()[0]

# Procedure: tree_clade_update_variants
# Purpose: replace clade's variants with a different set of variants
def tree_clade_update_variants(dbo, clade, variants):
    dc = dbo.cursor()
    dc.execute('delete from cladevariants where cladeid=?', (clade,))
    newlist = [(clade, v) for v in variants]
    dc.executemany('insert into cladevariants values(?,?)', (newlist))
    dbo.commit()
    return

# Procedure: tree_clade_update_kits
# Purpose: replace clade's kits with a different set of kits
def tree_clade_update_kits(dbo, clade, kits):
    dc = dbo.cursor()
    dc.execute('delete from cladekits where cladeid=?', (clade,))
    newlist = [(clade, v) for v in kits]
    dc.executemany('insert into cladekits values(?,?)', (newlist))
    dbo.commit()
    return

# Procedure: tree_uptree_variants
# Purpose: get all of the variants of the parent nodes
def tree_uptree_variants(dbo, clade):
    ancestors = tree_clade_get_ancestors(dbo, clade)
    vset = set()
    for a in ancestors:
        vset.update(set(tree_clade_get_variants(dbo, a)))
    return vset

# Procedure: tree_clade_variants_in_ancestors
# Purpose: check ancestors of a clade for a set of variants
def tree_clade_variants_in_ancestors(dbo, clade, variants):
    vset = tree_uptree_variants(dbo, clade)
    if set(variants).issubset(vset):
        return True
    return False

# Procedure: tree_add_child_clade
# Purpose: insert a new child of a specified node
def tree_add_child_clade(dbo, parent, newchild):
    trace(3, 'tree_add_child_clade adds {} under {}'.format(newchild,parent))
    trace(9, '  with variants {}'.format(tree_clade_get_variants(dbo,newchild)))
    # no work to do - already a child
    # this is useful if we want to call this more than once for the same child
    if newchild in tree_clade_get_children(dbo, parent):
        return
    dc = dbo.cursor()
    # all parents point to child node
    dc.execute('''INSERT INTO treepaths
                     (ancestor, descendant, treedepth)
                  SELECT ancestor, ?, treedepth+1 from treepaths
                  where descendant=?''', (newchild, parent))
    dbo.commit()
    # child clade contains variants of parent - remove them from child
    pvars = set(tree_clade_get_variants(dbo, parent))
    cvars = set(tree_clade_get_variants(dbo, newchild))
    newvars = cvars - pvars
    samevars = pvars.intersection(cvars)
    if not pvars:
        trace(3, 'no variants: parent node {}'.format(parent))
    if samevars:
        tree_clade_update_variants(dbo, parent, samevars)
    if newvars:
        tree_clade_update_variants(dbo, newchild, newvars)

    # child clade contains kits of parent - remove them from parent
    pkits = set(tree_clade_get_kits(dbo, parent))
    ckits = set(tree_clade_get_kits(dbo, newchild))
    tree_clade_update_kits(dbo, parent, pkits - ckits)
    dbo.commit()
    return

# Procedure: tree_clade_get_ancestors
# Purpose: return cladeIDs for all ancestors of a clade
# Return: list of ancestor IDs, including the child
def tree_clade_get_ancestors(dbo, child):
    dc = dbo.cursor()
    dc.execute('''SELECT c.id FROM treeclade c
                  JOIN treepaths t on c.ID=t.ancestor
                  WHERE t.descendant=? and t.treedepth > 0
                  ORDER BY treedepth''', (child,))
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
                  WHERE t.ancestor=? and t.treedepth=?''', (parent,1))
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
def tree_insert_clade(dbo, newclade, parent, child):
    dc = dbo.cursor()

    # remove old paths to child, add child under parent
    dc.execute('delete from treepaths where descendant=? and treedepth>0',
                   (newclade,))
    dbo.commit()
    tree_add_child_clade(dbo, parent, newclade)
    dbo.commit()
    # remove old paths to child, add child under parent
    dc.execute('delete from treepaths where descendant=? and treedepth>0',
                   (child,))
    dbo.commit()
    # add previous child below newclade
    tree_add_child_clade(dbo, newclade, child)
    dbo.commit()
    # all children of previous child are now depth+1
    dc.execute('''update treepaths set treedepth=treedepth+1
                  where ancestor=? and descendant in
                  (select descendant from treepaths
                     where ancestor=? and treedepth>0)
                    and treedepth>0''', (parent,child))

    # all children of previous child are also a descendant of newclade
    dc.execute('''insert into treepaths(ancestor,descendant,treedepth)
                  select ?,descendant,treedepth+1 from treepaths
                  where ancestor=? and treedepth>0''', (newclade, child))
    dbo.commit()
    return


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
    tree_clade_update_kits(dbo, node, allkits)
    return

# Procedure: tree_clade_get_variants
# Purpose: return list of variant ids associated with node
def tree_clade_get_variants(dbo, node, include_ancestors=0):
    dc = dbo.cursor()
    nodes = tree_clade_get_ancestors(dbo, node)[:include_ancestors] + [node]
    trace(9, 'tree_clade_get_variants nodes = {}'.format(nodes))
    vidlist = []
    for clade in nodes:
        dc.execute('select vid from cladevariants where cladeid=?', (clade,))
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
    if not tree_is_leaf_node(dbo, node):
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

# Procedure: tree_add_children
# Purpose: copy children of a given node under a new node
def tree_add_children(dbo, fromnode, intonode, remainder):
    dc = dbo.cursor()
    for child in tree_clade_get_children(dbo, fromnode):
        if child in remainder:
            continue
        trace(5, 'adding child {} under {}'.format(child, intonode))
        # need to copy the node and its sub-tree, not just the node
        remainder.append(child)
        tree_clone_tree(dbo, child, intonode)
        dbo.commit()

# Procedure: tree_combine_children
# Purpose: collapse identical children under a clade into one
def tree_combine_children(dbo, node):
    children = tree_clade_get_children(dbo, node)
    for ii,child in enumerate(children):
        vars1 = set(tree_clade_get_variants(dbo, child))
        kits1 = set(tree_clade_get_kits(dbo, child))
        for other in children[ii+1:]:
            vars2 = set(tree_clade_get_variants(dbo, other))
            kits2 = set(tree_clade_get_kits(dbo, other))
            # two child nodes share some variants
            if vars1.intersection(vars2):
                trace(8, 'combining {} into {}'.format(other, child))
                remainder = []
                tree_merge_into(dbo, other, child, remainder)
                # FIXME - does loop we're in work if we delete the node?
                # I think it's OK - kits and variants would be empty?
                tree_delete_tree(dbo, other)
                dbo.commit()
    return

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
        remainder.append(fromtree)
        dbo.commit()
        return

    # process all of the children of the "from" tree
    for child in tree_clade_get_children (dbo, fromtree):
        trace(3, 'tree_merge_into - merge child: {}'.format(child))
        if tree_is_leaf_node(dbo, child):
            tree_clade_merge_into(dbo, child, intotree, remainder)
        else:
            tree_merge_into(dbo, child, intotree, remainder)
        remainder.append(child)

    # after processing all of this node's children, process this node
    if not tree_is_leaf_node(dbo, fromtree):
        tree_clade_merge_into(dbo, fromtree, intotree, remainder)
        remainder.append(fromtree)
    dbo.commit()
    return

# Procedure: tree_clade_merge_into
# Purpose: merge a node into another tree
# Info:
#   WORK IN PROGRESS
#   There are many issues with this merging of two trees.
#   A clade is defined as a set of variants that a set of kits has. To merge it
#   into another tree, the variants in common with some node of the other tree
#   need to be factored out, then the node needs to be merged under that node
#   that had some of the variants. If no variants are found in common anywhere
#   in the tree, the new node needs to be (perhaps) parented under the top node
#   of the tree. When some variants are shared with the parent, the parent
#   needs to be split into two, with the new parent being all of the shared
#   variants. Trees are not perfect (e.g. variants that "should" be called are
#   not, and this causes merging problems, which are currently not dealt with.
def tree_clade_merge_into(dbo, fromnode, intotree, remainder):

    if fromnode in remainder:
        return

    trace(2, 'merge {} into {}'.format(fromnode, intotree))

    fromvars = set(tree_clade_get_variants(dbo, fromnode))

    # skip over treetop node or other invalid nodes
    if not fromvars:
        trace(2, 'empty node skipped: {},{}'.format(fromnode,intotree))
        remainder.append(fromnode)
        return

    fromkits = tree_clade_get_kits(dbo, fromnode)

    # traverse the target tree from the bottom up (leaf nodes first)
    for node in tree_clades_depth_first(dbo, intotree):
        intovars =set(tree_clade_get_variants(dbo, node, 0))

        # Cases: 1) set of variants in "from" node identically matches the set
        # of variants of some node in "to" tree; 2) "from" set is a subset of
        # some node's variants in "to" tree; 3) there are no shared variants
        # with any node in "to" tree. NB: the algorithm assumes that we already
        # handled all nodes below (descendants of) fromnode. If there are no
        # shared variants with any node, the node may still be merged later on,
        # when a higher node is merged.  Fromnode might also represent a span
        # of nodes in the "to" tree. E.g. if "to" tree has a>b>c and "from"
        # node has a,b,c, we should still merge the from node. This is handled
        # as case 4) "from" node has a superset of the variants in a node of
        # "to" tree. Then we can look up the tree at the ancestors and ensure
        # the variants are above and merge if yes. Finally, 5) the "from" node
        # might have some of the variants of the "to" node, and we still need
        # to merge and create a new node.

        # CASE 1: variants are identical to some node - copy kits and be done
        if intovars == fromvars:
            trace(8, 'CASE 1')
            tree_clade_add_kits(dbo, node, fromkits)
            # add all of the node's children and subtree to the "into" tree
            tree_add_children(dbo, fromnode, node, remainder)
            remainder.append(fromnode)
            #tree_combine_children(dbo, node)
            dbo.commit()
            return

        # CASE 4: variants are a superset of some node - split clade
        elif intovars and fromvars.issuperset(intovars):
            trace(8, 'CASE 4')
            trunkvars = tree_uptree_variants(dbo, node)
            # there are always some shared variants in this CASE
            shared = fromvars.intersection(trunkvars.union(intovars))
            # if there are unshared vars, they'll create a new clade
            unshared = fromvars - shared
            trace(10, 'fromvars: {}; trunkvars: {}; intovars: {}'.format(
                fromvars, trunkvars, intovars))
            trace(10, 'shared: {}; unshared: {}'.format(shared, unshared))
            if unshared:
                newclade = tree_newclade(dbo, fromkits, unshared)
                tree_add_child_clade(dbo, node, newclade)
                bottom = newclade
            else:
                tree_clade_add_kits(dbo, node, fromkits)
                bottom = node
            # add all of the node's children and subtree to the "into" tree
            remainder.append(fromnode)
            tree_add_children(dbo, fromnode, bottom, remainder)
            #tree_combine_children(dbo, bottom)
            dbo.commit()
            return

        # CASE 2: variants are a subset of some node
        # No child of this node contains any since we traverse depth first.
        # this means shared variants need to be pushed down to new node
        elif fromvars and fromvars.issubset(intovars):
            trace(8, 'CASE 2')
            shared = fromvars
            unshared = intovars - fromvars
            # shared variants move down to the new clade
            tree_clade_update_variants(dbo, node, unshared)
            # all existing kits move down to the new clade
            kits = tree_clade_get_kits(dbo, node) + fromkits
            tree_clade_update_kits(dbo, node, [])
            children = tree_clade_get_children(dbo, node)
            newclade = tree_newclade(dbo, kits, shared)
            for child in children:
                trace(2, 'insert {} between {} and {}'.format(
                    newclade, node, child))
                # "insert" newclade between parent and child
                tree_insert_clade(dbo, newclade, node, child)
            remainder.append(fromnode)
            dbo.commit()
            return

        # CASE 5: some shared variants - need to pull out a new clade
        elif fromvars.intersection(intovars):
            trace(8, 'CASE 5')
            trunkvars = tree_uptree_variants(dbo, node)
            shared = fromvars.intersection(trunkvars)
            unshared = fromvars - shared
            newclade = tree_newclade(dbo, fromkits, unshared)
            tree_add_child_clade(dbo, node, newclade)
            remainder.append(fromnode)
            dbo.commit()

    # CASE 3
    # reaching here means no shared variants found in the tree
    trace(8, 'CASE 3')
    # punt for now - most of these are handled by an ancestor being merged
    return
    treetop = get_treetops(dbo, intotree)
    tree_add_child_clade(dbo, treetop, newclade)
    dbo.commit()

# Procedure: tree_to_dot
# Purpose: turn a tree into DOT language
# Info:
#   Returns a graphviz string, viewable with xdot or similar.  Rendering with
#   raw nodes (dereference=False) is mainly useful for debugging
#   if there are many variants at a node, compact=True omits printing them
def tree_to_dot(dbo, top, dereference=True, compact=False):
    def nodestr(nodeid, kits, variants):
        return '{} [label="{}"];\n'.format(
            str(nodeid), '\\n'.join(['{'+str(nodeid)+'}'] +
                    [str(v) for v in variants]+['kit '+str(k) for k in kits]))
    def linkstr(parentid, childid):
        return '{} -> {};\n'.format(str(parentid), str(childid))
    nodelist = tree_clade_get_descendants(dbo, top)
    if dereference:
        allvars = set()
        allkits = set()
        for node in nodelist:
            allvars.update(set(tree_clade_get_variants(dbo, node)))
            allkits.update(set(tree_clade_get_kits(dbo, node)))
        vardefs = get_variant_defs(dbo, allvars)
        kitdefs = get_kit_ids(dbo, allkits)

    out = ''
    for node in nodelist:
        variants = tree_clade_get_variants(dbo,node)
        kits = tree_clade_get_kits(dbo,node)
        nvars = len(variants)
        if dereference:
            kitlabels = [kitdefs[k] for k in kits]
            varlabels = [vardefs[v][3] for v in variants]
        else:
            kitlabels = list(kits)
            varlabels = variants

        if compact and nvars > 5:
            varlabels = ['{} variants'.format(nvars)]

        out += nodestr(node, kitlabels, varlabels)


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

    # purge all tree data, for all trees
    if False:
        dc = dbo.cursor()
        dc.execute('delete from treepaths')
        dc.execute('delete from treeclade')
        dc.execute('delete from cladekits')
        dc.execute('delete from cladevariants')
        dbo.commit()
        sys.exit(0)

    # delete all treetop trees; this may leave orphaned entries
    if False:
        for tree in tree_get_treetops(dbo):
            trace(0, 'delete tree {}'.format(tree))
            tree_delete_tree(dbo, tree)
        dbo.commit()
        sys.exit(0)

    # merge all treetop trees into lowest-numbered one
    if True:
        tree_cleanup(dbo)
        trees = sorted(tree_get_treetops(dbo), reverse=True)
        # index of tree for target of merge in sorted array of trees
        totree = -1
        # maximum number of trees to merge into it
        numtrees = 5
        for tree in trees[totree-numtrees:totree]:
            remainder = []
            tree_merge_into(dbo,tree,trees[totree],remainder)
            trace(0, 'merged tree {} into {}'.format(tree,trees[totree]))
            nodelist = tree_clade_get_descendants(dbo, trees[totree])
            trace(0, 'combining children of all nodes')
            for n in nodelist:
                tree_combine_children(dbo, n)
                dbo.commit()
            trace(0, 'merging identical nodes across the tree')
            nodelist = tree_clades_depth_first(dbo, trees[totree])
            nodes = enumerate(nodelist)
            # if variants are identical merge nodes - into the deepest node
            variants = [set(tree_clade_get_variants(dbo,n)) for n in nodelist]
            depth = [tree_node_depth(dbo,n) for n in nodelist]
            trace(0, 'variants: {}'.format(variants))
            trace(0, 'depth: {}'.format(depth))
            merge_into = []
            ndict = dict()
            for ii,node in nodes:
                ndict[ii] = node
                trace(0, 'node {} (i={})'.format(node, ii))
                for jj in range(ii+1, len(nodelist)):
                    if variants[ii] == variants[jj]:
                        if depth[ii] < 4 or depth[jj] < 4:
                            continue
                        if depth[ii] > depth[jj]:
                            merge_into.append((jj, ii))
                        else:
                            merge_into.append((ii, jj))
            trace(0, 'length of merge_into: {}'.format(len(merge_into)))

            dbo.commit()

            done = []
            trace(0, 'ndict: {}'.format(ndict))
            for ifrom,iinto in merge_into:
                if (ndict[ifrom] in done) or (ndict[iinto] in done):
                    trace(0, 'unable to merge {} into {}'.format(
                        ndict[ifrom], ndict[iinto]))
                    continue
                trace(0, 'cleanup node {} by merging into {}'.format(
                    ndict[ifrom], ndict[iinto]))
                remainder = []
                tree_merge_into(dbo, ndict[ifrom], ndict[iinto], remainder)
                for nid in tree_clade_get_descendants(dbo, ndict[ifrom]):
                    done.append(nid)
                done.append(ndict[ifrom])
                tree_delete_tree(dbo, ndict[ifrom])
                dbo.commit()

        trace(0, 'writing tree to tree.dot...')
        # render with raw node and kit IDs
        open('tree-raw.gv','w').write(tree_to_dot(dbo, trees[totree], False))
        # render with actual kit names and variant names
        open('tree.gv','w').write(tree_to_dot(dbo, trees[totree]))
        sys.exit(0)

    # various unit tests
    if True:
        # create a simple tree by adding and inserting nodes
        treetop1 = tree_newclade(dbo, [], [], 'Top')
        c1 = tree_newclade(dbo, [1], [1])
        tree_add_child_clade(dbo, treetop1, c1)
        c2 = tree_newclade(dbo, [2], [2])
        tree_add_child_clade(dbo, c1, c2)
        c3 = tree_newclade(dbo, [3], [3])
        tree_add_child_clade(dbo, c2, c3)
        c5 = tree_newclade(dbo, [5], [5])
        tree_add_child_clade(dbo, c3, c5)
        open('tree1a.dot','w').write(tree_to_dot(dbo, treetop1, False))
        c4 = tree_newclade(dbo, [4], [4])
        tree_insert_clade(dbo, c4, c3, c5)
        open('tree1b.dot','w').write(tree_to_dot(dbo, treetop1, False))

        # create a tree that we will merge into
        treetop2 = tree_newclade(dbo, [], [], 'Top')
        trace(0, 'treetop2: {}'.format(treetop2))
        c1 = tree_newclade(dbo, [10,20,30,40], [5,6,7,8])
        tree_add_child_clade(dbo, treetop2, c1)
        c2 = tree_newclade(dbo, [10,20], [9])
        tree_add_child_clade(dbo, c1, c2)
        open('tree2a.dot','w').write(tree_to_dot(dbo, treetop2, False))

        # create a tree that will be merged in
        treetop3 = tree_newclade(dbo, [], [], 'Top')
        trace(0, 'treetop3: {}'.format(treetop3))
        c1 = tree_newclade(dbo, [30,31,32], [5,6,7,8])
        tree_add_child_clade(dbo, treetop3, c1)
        c2 = tree_newclade(dbo, [30,31], [9])
        tree_add_child_clade(dbo, c1, c2)
        c4 = tree_newclade(dbo, [50], [10])
        tree_add_child_clade(dbo, c2, c4)
        c3 = tree_newclade(dbo, [10,20,40,41], [5,6,7])
        tree_insert_clade(dbo, c3, treetop3, c1)
        open('tree3a.dot','w').write(tree_to_dot(dbo, treetop3, False))

        # merge tree 3 into tree 2
        remainder = []
        tree_merge_into(dbo,treetop3,treetop2,remainder)

        open('tree2b.dot','w').write(tree_to_dot(dbo, treetop2, False))

        # call print on each node of tree 1
        tree_func_bottom_up(dbo, treetop1, print)
        trace(0, '{}'.format(tree_clades_depth_first(dbo,treetop1)))
        trace(0, '{}'.format(tree_clade_get_ancestors(dbo,treetop1)))

        tree_delete_tree(dbo, treetop1)
        tree_delete_tree(dbo, treetop2)
        tree_delete_tree(dbo, treetop3)

        sys.exit(0)
