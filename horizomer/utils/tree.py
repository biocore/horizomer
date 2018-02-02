#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

#
# functions relevant to phylogenetic tree operations
#

from skbio import TreeNode


def support(node):
    """Get support value of a node.

    Parameters
    ----------
    node : skbio.TreeNode
        node to get support value of

    Returns
    -------
    float or None
        support value of the node, or None if not available

    Notes
    -----
    A "support value" is defined as the numeric form of a whole node label
    without ":", or the part preceding the first ":" in the node label.
    - For examples: "(a,b)1.0", "(a,b)1.0:2.5", and "(a,b)'1.0:species_A'". In
    these cases the support values are all 1.0.
    - For examples: "(a,b):1.0" and "(a,b)species_A". In these cases there are
    no support values.
    """
    try:
        return float(node.name.split(':')[0])
    except (ValueError, AttributeError):
        return None


def unpack(node):
    """Unpack an internal node of a tree.

    Parameters
    ----------
    node : skbio.TreeNode
        node to unpack

    Notes
    -----
    This function sequentially: 1) elongates child nodes by branch length of
    self (omit if there is no branch length), 2) removes self from parent node,
    and 3) grafts child nodes to parent node.

    Here is an illustration of the "unpack" operation:
                /----a
          /c---|
         |      \--b
    -----|
         |        /---d
          \f-----|
                  \-e

    Unpack node "c" and the tree becomes:
          /---------a
         |
    -----|--------b
         |
         |        /---d
          \f-----|
                  \-e
    """
    if node.is_root():
        raise ValueError('Cannot unpack root.')
    parent = node.parent
    blen = (node.length or 0.0)
    for child in node.children:
        clen = (child.length or 0.0)
        child.length = (clen + blen or None)
    parent.remove(node)
    parent.extend(node.children)


def has_duplicates(tree):
    """Test whether there are duplicated taxa (tip names) in a tree.

    Parameters
    ----------
    tree : skbio.TreeNode
        tree for check for duplicates

    Returns
    -------
    bool
        whether there are duplicates
    """
    taxa = [tip.name for tip in tree.tips()]
    if '' in taxa or None in taxa:
        raise ValueError('Empty taxon name(s) found.')
    return len(set(taxa)) < len(taxa)


def compare_topology(tree1, tree2):
    """Test whether the topologies of two trees with all nodes assigned
    unique IDs are identical.

    Parameters
    ----------
    tree1 : skbio.TreeNode
        first tree in comparison
    tree2 : skbio.TreeNode
        second tree in comparison

    Returns
    -------
    bool
        whether the topologies are identical

    Notes
    -----
    Given all nodes (internal nodes or tips) have unique IDs, one just needs
    to check if all node-to-parent pairs are identical.
    Flipping child nodes does not change topology. Branch lengths are ignored
    when comparing topology.
    """
    n2p1, n2p2 = ({node.name: node.parent.name
                  for node in tree.traverse() if not node.is_root()}
                  for tree in (tree1, tree2))
    return n2p1 == n2p2


def intersect_trees(tree1, tree2):
    """Shrink two trees to contain only overlapping taxa.

    Parameters
    ----------
    tree1 : skbio.TreeNode
        first tree to intersect
    tree2 : skbio.TreeNode
        second tree to intersect

    Returns
    -------
    tuple of two TreeNodes
        resulting trees containing only overlapping taxa
    """
    for tree in (tree1, tree2):
        if has_duplicates(tree):
            raise ValueError('Either tree has duplicated taxa.')
    taxa1 = set([tip.name for tip in tree1.tips()])
    taxa2 = set([tip.name for tip in tree2.tips()])
    taxa_lap = taxa1.intersection(taxa2)
    if len(taxa_lap) == 0:
        raise ValueError('Trees have no overlapping taxa.')
    tree1_lap = tree1.shear(taxa_lap)
    tree2_lap = tree2.shear(taxa_lap)
    return (tree1_lap, tree2_lap)


def unpack_short_branch_nodes(tree, cutoff):
    """Unpack internal nodes with branch length below cutoff.

    Parameters
    ----------
    tree : skbio.TreeNode
        tree to search for nodes to unpack
    cutoff : int or float
        branch length cutoff under which node should be unpackd

    Returns
    -------
    skbio.TreeNode
        resulting tree with short branches unpackd

    Notes
    -------
    Nodes without branch length are considered as having zero branch length and
    will be unpackd.
    """
    tcopy = tree.copy()
    nodes_to_unpack = []
    for node in tcopy.non_tips():
        if (node.length or 0.0) < cutoff:
            nodes_to_unpack.append(node)
    for node in nodes_to_unpack:
        unpack(node)
    return tcopy


def unpack_low_support_nodes(tree, cutoff):
    """Unpack internal nodes with support value below cutoff.

    Parameters
    ----------
    tree : skbio.TreeNode
        tree to search for nodes to unpack
    cutoff : int or float
        node support value cutoff under which it should be unpackd

    Returns
    -------
    skbio.TreeNode
        resulting tree with low-support nodes unpackd

    Notes
    -------
    Nodes without support value will NOT be unpackd.
    """
    tcopy = tree.copy()
    nodes_to_unpack = []
    for node in tcopy.non_tips():
        spt = support(node)
        if spt is not None and spt < cutoff:
            nodes_to_unpack.append(node)
    for node in nodes_to_unpack:
        unpack(node)
    return tcopy


def read_taxdump(nodes_fp, names_fp=None):
    """Read NCBI taxdump.

    Parameters
    ----------
    nodes_fp : str
        file path to NCBI nodes.dmp
    names_fp : str, optional
        file path to NCBI names.dmp

    Returns
    -------
    dict of dict
        taxid : {
            'parent' : str
                parent taxid
            'rank' : str
                taxonomic rank
            'name' : str
                taxon name, empty if names_fp is None
            'children' : set of str
                child taxids
        }
    """
    taxdump = {}

    # format of nodes.dmp: taxid | parent taxid | rank | more info...
    with open(nodes_fp, 'r') as f:
        for line in f:
            x = line.rstrip('\r\n').replace('\t|', '').split('\t')
            taxdump[x[0]] = {'parent': x[1], 'rank': x[2], 'name': '',
                             'children': set()}

    # format of names.dmp: taxid | name | unique name | name class |
    if names_fp is not None:
        with open(names_fp, 'r') as f:
            for line in f:
                x = line.rstrip('\r\n').replace('\t|', '').split('\t')
                if x[3] == 'scientific name':
                    taxdump[x[0]]['name'] = x[1]

    # identify children taxids
    for tid in taxdump:
        pid = taxdump[tid]['parent']
        if tid != pid:  # skip root whose taxid equals its parent
            taxdump[pid]['children'].add(tid)

    return taxdump


def build_taxdump_tree(taxdump):
    """Build NCBI taxdump tree.

    Parameters
    ----------
    taxdump : dict of dict
        attributes of each taxid, see read_taxdump

    Returns
    -------
    skbio.TreeNode
        a tree representing taxdump
    """
    # create the tree from root
    tree = TreeNode('1')

    # iteratively attach child nodes to parent node
    def iter_node(node):
        for cid in taxdump[node.name]['children']:
            child = TreeNode(cid)
            node.extend([child])
            iter_node(child)

    iter_node(tree)
    return tree
