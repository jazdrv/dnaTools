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
from functools import reduce

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
    dc.execute('''delete from treepaths where descendant not in
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
def tree_clone_tree(dbo, fromtree, destination, remainder):
    nvars = tree_clade_get_variants(dbo, fromtree)
    nkits = tree_clade_get_kits(dbo, fromtree)
    newchild = tree_newclade(dbo, nkits, nvars)
    remainder.update(nvars)
    tree_add_child_clade(dbo, destination, newchild)
    dbo.commit()
    for child in tree_clade_get_children(dbo, fromtree):
        tree_clone_tree(dbo, child, newchild, remainder)
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
    try:
        dc.execute('''INSERT INTO treepaths (ancestor, descendant, treedepth)
                      values(?,?,0)''', (nodeID, nodeID))
    except:
        trace(0, 'path insert failed for node {}'.format(nodeID))
        raise

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
    # all parents point to new child node and its descendants
    trace(8, 'tree link maintenance, parent={}'.format(parent))
    dc.execute('select * from treepaths')
    for tup in dc:
        trace(20, 'treepath {}'.format(tup))
    dc.execute('''INSERT or ignore INTO treepaths
                        (ancestor, descendant, treedepth)
                  select a.ancestor, b.descendant, a.treedepth+b.treedepth+1
                  from treepaths a, treepaths b
                  where a.descendant=? and b.ancestor=?''', (parent,newchild))
    dbo.commit()
    trace(8, 'after...'.format(parent))
    dc.execute('select * from treepaths')
    for tup in dc:
        trace(20, 'treepath {}'.format(tup))
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

    # remove old paths to child
    dc.execute('''delete from treepaths
                    where ancestor in
         (select ancestor from treepaths where treedepth>0 and descendant=?)
                    and descendant in
         (select descendant from treepaths where treedepth>0 and ancestor=?)
                   ''', (child, parent))
    dbo.commit()
    # add child below new clade
    tree_add_child_clade(dbo, newclade, child)
    # add new clade below parent
    tree_add_child_clade(dbo, parent, newclade)
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
    allkits = set(list(nodekits) + list(newkits))
    tree_clade_update_kits(dbo, node, allkits)
    return

# Procedure: tree_clade_add_variants
# Purpose: add a list of variants to a given node
def tree_clade_add_variants(dbo, node, newvariants):
    nodevariants = tree_clade_get_variants(dbo, node)
    allvariants = set(list(nodevariants) + list(newvariants))
    tree_clade_update_variants(dbo, node, allvariants)
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

# Procedure: tree_get_leaf_nodes
# Purpose: find all leaf nodes in a tree
def tree_get_leaf_nodes(dbo, node):
    if tree_is_leaf_node(dbo, node):
        return [node]
    else:
        leafs = []
        for child in tree_clade_get_children(dbo, node):
            leafs += tree_get_leaf_nodes(dbo, child)
    return leafs

# Procedure: tree_get_paths
# Purpose: find all paths from leaf nodes to top of tree
def tree_get_paths(dbo, node):
    leafnodes = tree_get_leaf_nodes(dbo, node)
    paths = []
    for leafnode in leafnodes:
        variants = []
        kits = []
        nodelist = [leafnode] + tree_clade_get_ancestors(dbo, leafnode)
        for treenode in nodelist:
            variants.append(set(tree_clade_get_variants(dbo, treenode)))
            kits.append(set(tree_clade_get_kits(dbo, treenode)))
        paths.append((nodelist, kits, variants))
    return paths

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
        cvars = set(tree_clade_get_variants(dbo, child))
        if cvars.issubset(remainder):
            continue
        trace(5, 'adding child {} under {}'.format(child, intonode))
        # need to copy the node and its sub-tree, not just the node
        remainder.update(cvars)
        tree_clone_tree(dbo, child, intonode, remainder)
        dbo.commit()

# Procedure: tree_combine_children
# Purpose: collapse identical children under a clade into one
def tree_combine_children(dbo, node):
    children = tree_clade_get_children(dbo, node)
    visited = []
    for ii,child in enumerate(children):
        if child in visited:
            continue
        visited.append(child)
        vars1 = set(tree_clade_get_variants(dbo, child))
        if not vars1:
            # FIXME should not happen
            continue
        for other in children[ii+1:]:
            if other in visited:
                continue
            vars2 = set(tree_clade_get_variants(dbo, other))
            remainder = set()
            # two child nodes share some variants
            if vars1 == vars2:
                tree_clade_add_kits(dbo, child, tree_clade_get_kits(dbo, other))
                for sub in tree_clade_get_children(dbo, other):
                    tree_clone_tree(dbo, sub, child, remainder)
                tree_delete_tree(dbo, other)
                visited.append(other)
                dbo.commit()
            elif vars1.issubset(vars2):
                trace(8, 'combining {} into {}'.format(other, child))
                tree_clade_update_variants(dbo, other, vars2-vars1)
                tree_clone_tree(dbo, other, child, remainder)
                tree_delete_tree(dbo, other)
                visited.append(other)
                dbo.commit()
            elif vars2.issubset(vars1):
                kits1 = tree_clade_get_kits(dbo, child)
                kits2 = tree_clade_get_kits(dbo, other)
                tree_clade_update_kits(dbo, child, kits2)
                tree_clade_update_variants(dbo, child, vars2)
                newnode = tree_newclade(dbo, kits1, vars1-vars2)
                tree_add_child_clade(dbo, child, newnode)
                for sub in tree_clade_get_children(dbo, other):
                    tree_clone_tree(dbo, sub, child, remainder)
                tree_delete_tree(dbo, other)
                visited.append(other)
                vars1 = vars2
                dbo.commit()
    return

# Procedure: path_compare
# Purpose: see if there is sufficient variant overlap to merge tree nodes
def path_compare(fromsets, intosets):
    if not fromsets[0] or not intosets[0]:
        return False
    if not fromsets[0].intersection(intosets[0]):
        return False
    fvars = set()
    tvars = set()
    [fvars.update(s) for s in fromsets]
    [tvars.update(s) for s in intosets]
    overlap = fvars.intersection(tvars)
    if len(overlap) < 2:
        return False
    # FIXME - maybe want some path checking, more rigorous check?
    return True


# Procedure: tree_combine_orphans
# Purpose: eliminate empty nodes and combine identical nodes
def tree_combine_orphans(dbo, tree):
    nodes = tree_clades_depth_first(dbo, tree)
    for node in nodes:
        remainder = set()
        nvars = set(tree_clade_get_variants(dbo, node))
        nkits = tree_clade_get_kits(dbo, node)
        try:
            pnode = tree_clade_get_parent(dbo, node)
        except:
            trace(0, 'did not fetch parent of {}'.format(node))
            continue
        if (not nvars) or (nvars == set(tree_clade_get_variants(dbo, pnode))):
            for child in tree_clade_get_children(dbo, node):
                tree_clone_tree(dbo, child, pnode, remainder)
            tree_delete_tree(dbo, node)
            tree_clade_add_kits(dbo, pnode, nkits)
            dbo.commit()
            
# Procedure: tree_combine_sibs
# Purpose: eliminate duplicate parallel nodes
def tree_combine_sibs(dbo, tree):
    nodes = tree_clades_depth_first(dbo, tree)
    for node in nodes:
        if tree_is_leaf_node(dbo, node):
            continue
        children = tree_clade_get_children(dbo, node)
        if len(children) > 1:
            tree_combine_children(dbo, node)
            dbo.commit()
    return

# Procedure: tree_pushup_variants
# Purpose: move variants higher on tree where shared on branches
def tree_pushup_variants(dbo, tree):
    paths = tree_get_paths(dbo, tree)
    dc = dbo.cursor()
    dc.execute('select vid, cladeid from cladevariants order by 1')
    vid1,clade1 = dc.fetchone()
    for vid2, clade2 in dc:
        if vid2 == vid1:
            path1 = [p[0] for p in paths if clade1 in p[0]]
            path2 = [p[0] for p in paths if clade2 in p[0]]
            if not path1 or not path2:
                trace(0, 'WARN in pushup {}.{}, {}.{}'.format(
                    clade1,vid1, clade2,vid2))
                continue
            # push the variant up to the highest shared node
            pshared = enumerate(zip(reversed(path1[0]), reversed(path2[0])))
            for ii, (n1,n2) in pshared:
                if n1 == n2:
                    shared_node = n1
                else:
                    break
            if ii > 0:
                vidset = set([vid1])
                tree_clade_add_variants(dbo, shared_node, vidset)
                newvars = set(tree_clade_get_variants(dbo, clade1)) - vidset
                tree_clade_update_variants(dbo, clade1, newvars)
                newvars = set(tree_clade_get_variants(dbo, clade2)) - vidset
                tree_clade_update_variants(dbo, clade2, newvars)
        else:
            vid1=vid2
            clade1=clade2
    return

# Procedure: tree_merge_into
# Purpose: merge one tree into another tree
# Info:
#   The "into" tree should be the top node of a tree to receive the new tree
#   data. The "from" tree should be a node that we will attempt to merge into
#   the receiving tree. When a node can't be merged, it's added to the
#   remainder tree.
#@profile
def tree_merge_into(dbo, fromtree, intotree, remainder):

    ipaths = None
    fpaths = tree_get_paths(dbo, fromtree)

    unhandled = set()
    # path: (nodes, kits, variants)
    for ii,path in enumerate(fpaths):

        if not ipaths:
            ipaths = tree_get_paths(dbo, intotree)

        trace(8, 'from-tree: paths: {}'.format(fpaths))
        trace(8, 'into-tree: paths: {}'.format(ipaths))

        # list of all vars on this path
        pathvars = reduce(lambda x,y: x.union(y), path[2])
        best_overlap = 0
        # count of variants in common with the path on the "from" tree
        path_overlap = [len(pathvars.intersection(s)) for s in
                            [reduce(lambda x,y: x.union(y), inp[2])
                                 for inp in ipaths]]
        trace(8, 'path-overlap: {}'.format(path_overlap))

        if path_overlap:
            best_overlap = max(path_overlap)
        trace(8, 'best path count {}: {}'.format(best_overlap,
                                ipaths[path_overlap.index(best_overlap)]))

        #if best_overlap == 0:
        if best_overlap < 3:
            # add under top of tree
            unhandled = unhandled - set(path[0])
            tree_clone_tree(dbo, path[0][-1], intotree, remainder)
            ipaths = None

        else:

            ipath = ipaths[path_overlap.index(best_overlap)]

            unhandled = unhandled - set(ipath[0])

            tree_clade_merge_into(dbo, path, ipath, remainder)
            ipaths = None

    if unhandled:
        trace(0, 'UNHANDLED NODES: {}'.format(unhandled))

    dbo.commit()
    tree_cleanup(dbo)
    #tree_pushup_variants(dbo, intotree)
    tree_combine_orphans(dbo, intotree)
    tree_combine_sibs(dbo, intotree)

    return


def tree_get_nodes(dbo, treetop):
    nodes = []
    children = tree_clade_get_children(dbo, treetop)
    for child in children:
        nodes += tree_get_nodes(dbo, child)
    nodes += [treetop]
    return nodes

# perform some sanity checking on tree
def sanity_check(dbo, tree):
    sane = True
    nodes = tree_get_nodes(dbo, tree)
    for ii, node in enumerate(nodes):
        for node2 in nodes[ii+1:]:
            s1 = tree_clade_get_variants(dbo, node)
            s2 = tree_clade_get_variants(dbo, node2)
            if s1 == s2:
                trace(1, 'identical nodes {},{} with variants {},{}'.format(
                    node, node2, s1, s2))
                sane = False
            if not s1:
                trace(1, 'empty node {}'.format(node))
                sane = False
    return sane

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
def tree_clade_merge_into(dbo, frompath, intopath, remainder):

    # path[0] is the leaf node; path[-1] is the top
    fromnodes = frompath[0]
    fromkits = frompath[1]
    fromvars = frompath[2]
    intonodes = intopath[0]
    intokits = intopath[1]
    intovars = intopath[2]

    trace(3, 'merge {} into {}'.format(frompath[0], intopath[0]))
    trace(8, 'fromvars: {}, intovars: {}'.format(fromvars,intovars))
    trace(8, 'fromkits: {}, intokits: {}'.format(fromkits,intokits))

    # sanity check - all input nodes should have var(s) except perhaps topnode
    if len(fromvars) > 1 and not min(map(len, fromvars[:-1])):
        trace(2, 'empty node(s) skipped: {},{}'.format(frompath,intopath))
        return None

    # set of variants that have already been processed
    remainder = set()

    # traverse the target tree from the bottom up (leaf nodes first)
    for ff,fnode in enumerate(fromnodes):

        # for debugging
        if False and not sanity_check(dbo, intonodes[-1]):
            trace(1, 'SANITY check failed')

        # treetop might not have any variants - done
        if not fromvars[ff]:
            break

        intopath = tree_get_paths(dbo, intopath[0][0])[0]
        intonodes = intopath[0]
        intokits = intopath[1]
        intovars = intopath[2]
        for ii,inode in enumerate(intonodes):

            # Cases: 1) set of variants in "from" node identically matches the
            # set of variants of some node in "to" tree; 2) "from" set is a
            # subset of some node's variants in "to" tree; 3) there are no
            # shared variants with any node in "to" tree. NB: the algorithm
            # assumes that we already handled all nodes below (descendants of)
            # fromnode. If there are no shared variants with any node, the node
            # may still be merged later on, when a higher node is merged.
            # Fromnode might also represent a span of nodes in the "to"
            # tree. E.g. if "to" tree has a>b>c and "from" node has a,b,c, we
            # should still merge the from node. This is handled as case 4)
            # "from" node has a superset of the variants in a node of "to"
            # tree. Then we can look up the tree at the ancestors and ensure
            # the variants are above and merge if yes. Finally, 5) the "from"
            # node might have some of the variants of the "to" node, and we
            # still need to merge and create a new node.

            # CASE 1: variants are identical to some node - copy kits, done
            if intovars[ii] == fromvars[ff]:
                trace(8, 'CASE 1')
                tree_clade_add_kits(dbo, inode, fromkits[ff])
                # add all of the node's children and subtree to the "into" tree
                tree_add_children(dbo, fnode, inode, remainder)
                tree_combine_children(dbo, inode)
                remainder.update(fromvars[ff])
                if ff:
                    remainder.update(fromvars[ff-1])
                dbo.commit()

            # CASE 4: variants are a superset of some node - add new clade
            # under the node. Original node keeps shared variants and
            # children. New node gets the unshared variants and the children of
            # fromnode.
            elif intovars[ii] and fromvars[ff].issuperset(intovars[ii]):
                trace(8, 'CASE 4')
                trunkvars = tree_uptree_variants(dbo, inode)
                # there are always some shared variants in this CASE
                shared = fromvars[ff].intersection(trunkvars.union(intovars[ii]))
                # if there are unshared vars, they'll create a new clade
                unshared = fromvars[ff] - shared
                trace(10, 'fromvars: {}; trunkvars: {}; intovars: {}'.format(
                    fromvars[ff], trunkvars, intovars[ii]))
                if unshared:
                    kits = tree_clade_get_kits(dbo, fnode)
                    newclade = tree_newclade(dbo, kits, unshared)
                    tree_add_child_clade(dbo, inode, newclade)
                    tree_add_children(dbo, fnode, newclade, remainder)
                    tree_combine_children(dbo, inode)
                    remainder.update(fromvars[ff])
                    if ff:
                        remainder.update(fromvars[ff-1])
                # it's possible all "new" variants are up-tree
                else:
                    tree_clade_add_kits(dbo, inode, fromkits[ff])
                    tree_add_children(dbo, fnode, inode, remainder)
                    tree_combine_children(dbo, inode)
                dbo.commit()

            # CASE 2: variants are a subset of some node. No child of this node
            # contains any overlap since we traverse depth first. A new node is
            # inserted between the node and its children and the shared
            # variants need to stay with original node. Unshared variants and
            # children go with the new node.
            elif fromvars[ff] and fromvars[ff].issubset(intovars[ii]):
                trace(8, 'CASE 2')
                trace(8, 'fromvars[ff]: {}, intovars[ii]: {}'.format(
                    fromvars[ff], intovars[ii]))
                trace(8, 'inode: {}, fnode: {}'.format(inode,fnode))
                trunkvars = tree_uptree_variants(dbo, fnode)
                # shared = variants in entire "from" tree also in "into" node
                shared = intovars[ii].intersection(trunkvars.union(fromvars[ff]))
                # subset: "from" node splits up the "into" node
                subvars = fromvars[ff].intersection(intovars[ii])
                # unshared means "into" node has variants not in "from" tree
                # the position of unshared vars is ambigous on the merged tree
                unshared = intovars[ii] - shared
                # up the tree: variants in both trees but higher up
                upvars = intovars[ii] - subvars - unshared
                trace(3, 'shared:{}'.format(shared))
                trace(3, 'unshared:{}'.format(unshared))
                trace(3, 'subvars:{}'.format(subvars))
                trace(3, 'upvars:{}'.format(upvars))


                # node at the top, variants in both trees
                trace(3, 'update node {} variants with {}'.format(inode, upvars))
                tree_clade_update_variants(dbo, inode, upvars)
                tree_clade_update_kits(dbo, inode, [])

                # variants lower in "from" tree than upvars
                newclade = tree_newclade(dbo, fromkits[ff], subvars)

                # FIXME - heuristics to place these unknowns
                # for now, place below subvars
                if unshared:
                    unshclade = tree_newclade(dbo, intokits[ii], unshared)
                else:
                    unshclade = None
                children = tree_clade_get_children(dbo, inode)
                for child in children:
                    trace(3, 'insert {} between {} and {}'.format(
                        newclade, inode, child))
                    # insert between parent and child of original node
                    tree_insert_clade(dbo, newclade, inode, child)
                    if unshclade:
                        tree_insert_clade(dbo, unshclade, newclade, child)

                remainder.update(subvars)

                # children of new node go under newly created clade
                tree_add_children(dbo, fnode, newclade, remainder)
                tree_combine_children(dbo, newclade)

                dbo.commit()

            # CASE 5: some shared variants - need to pull out a new clade
            elif fromvars[ff].intersection(intovars[ii]):
                trace(8, 'CASE 5')
                trunkvars = tree_uptree_variants(dbo, inode)
                shared = fromvars[ff].intersection(trunkvars)
                unshared1 = fromvars[ff] - shared
                tree_clade_update_variants(dbo, inode, shared)
                newclade1 = tree_newclade(dbo, fromkits[ff], unshared1)
                tree_add_child_clade(dbo, inode, newclade1)
                tree_add_children(dbo, fnode, newclade1, remainder)
                # new node for unshared variants of original node
                kits = tree_clade_get_kits(dbo, inode)
                tree_clade_update_kits(dbo, inode, [])
                unshared2 = intovars[ii] - shared
                newclade2 = tree_newclade(dbo, kits, unshared2)
                children = tree_clade_get_children(dbo, inode)
                for child in children:
                    trace(2, 'insert {} between {} and {}'.format(
                        newclade2, inode, child))
                    # "insert" newclade between parent and child of original node
                    tree_insert_clade(dbo, newclade2, inode, child)
                remainder.update(fromvars[ff])
                dbo.commit()

    # all nodes of frompath should be handled
    if False: #FIXME
        trace(0, 'WARN: not all nodes merged')
    return None

# Procedure: tree_to_dot
# Purpose: turn a tree into DOT language
# Info:
#   Returns a graphviz string, viewable with xdot or similar.  Rendering with
#   raw nodes (dereference=False) is mainly useful for debugging
#   if there are many variants at a node, compact=True omits printing them
def tree_to_dot(dbo, top, dereference=True, compact=False):
    def nodestr(nodeid, kits, variants):
        fontcolor=''
        if variants and type(variants[0]) == type('') and 'variants' in variants[0] \
                    and int(variants[0].split()[0]) > 30:
            fontcolor = 'fontcolor=blue'
        return '{} [{} label="{}"];\n'.format(
            str(nodeid), fontcolor, '\\n'.join(['{'+str(nodeid)+'}'] +
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

        if compact and nvars > 12:
            varlabels = ['{} variants'.format(nvars)]
        if compact and len(kits) > 12:
            kitlabels = ['{} kits'.format(len(kits))]

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

    if False:
        trace(0, 'leaf nodes: {}'.format(tree_get_leaf_nodes(dbo, 1)))
        sys.exit(0)

    # purge all tree data, for all trees
    if True:
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
    if False:
        tree_cleanup(dbo)
        trees = sorted(tree_get_treetops(dbo), reverse=True)
        trees = [5396, 5243, 5099]
        trees = [5533, 5396, 5243]
        trees = [5817, 5671, 5533, 5396]
        trees = [6099, 5960, 5817, 5671, 5533]
        trees = [6245, 6099, 5960, 5817, 5671]
        trees = [6529, 6385, 6245, 6099, 5960, 5817]
        trees = [6811, 6666, 6529, 6385, 6245, 6099, 5960]
        # index of tree for target of merge in sorted array of trees
        totree = -1
        # maximum number of trees to merge into it
        numtrees = 6
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
        # render compact with actual kit names and variant names
        open('tree.gv','w').write(tree_to_dot(dbo, trees[totree], True, True))
        # render with raw node and kit IDs
        open('tree-raw.gv','w').write(tree_to_dot(dbo, trees[totree], False))
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
        trace(0, 'leaf node: {}'.format(tree_get_leaf_nodes(dbo, treetop1)))
        trace(0, 'paths: {}'.format(tree_get_paths(dbo, treetop1)))
        open('tree1a.dot','w').write(tree_to_dot(dbo, treetop1, False))
        c4 = tree_newclade(dbo, [4], [4])
        tree_insert_clade(dbo, c4, c3, c5)
        trace(0, 'leaf node: {}'.format(tree_get_leaf_nodes(dbo, treetop1)))
        trace(0, 'paths: {}'.format(tree_get_paths(dbo, treetop1)))
        open('tree1b.dot','w').write(tree_to_dot(dbo, treetop1, False))

        # create a tree that we will merge into
        treetop2 = tree_newclade(dbo, [], [], 'Top')
        trace(0, 'treetop2: {}'.format(treetop2))
        c1 = tree_newclade(dbo, [10,20,30,40], [5,6,7,8])
        tree_add_child_clade(dbo, treetop2, c1)
        c2 = tree_newclade(dbo, [10,20], [9])
        tree_add_child_clade(dbo, c1, c2)
        trace(0, 'leaf node: {}'.format(tree_get_leaf_nodes(dbo, treetop2)))
        trace(0, 'paths: {}'.format(tree_get_paths(dbo, treetop2)))
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
        trace(0, 'leaf node: {}'.format(tree_get_leaf_nodes(dbo, treetop3)))
        trace(0, 'paths: {}'.format(tree_get_paths(dbo, treetop3)))
        open('tree3a.dot','w').write(tree_to_dot(dbo, treetop3, False))

        # sys.exit(0)

        # merge tree 3 into tree 2
        remainder = []
        tree_merge_into(dbo,treetop3,treetop2,remainder)

        trace(0, 'leaf node: {}'.format(tree_get_leaf_nodes(dbo, treetop2)))
        trace(0, 'paths: {}'.format(tree_get_paths(dbo, treetop2)))
        open('tree2b.dot','w').write(tree_to_dot(dbo, treetop2, False))

        # call print on each node of tree 1
        tree_func_bottom_up(dbo, treetop1, print)
        trace(0, '{}'.format(tree_clades_depth_first(dbo,treetop1)))
        trace(0, '{}'.format(tree_clade_get_ancestors(dbo,treetop1)))

        tree_delete_tree(dbo, treetop1)
        tree_delete_tree(dbo, treetop2)
        tree_delete_tree(dbo, treetop3)

        sys.exit(0)
