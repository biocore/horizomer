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
    taxa1 = set([tip.name for tip in tree1.tips()])
    taxa2 = set([tip.name for tip in tree2.tips()])
    taxa_lap = taxa1.intersection(taxa2)
    if len(taxa_lap) == 0:
        raise ValueError('Trees have no overlapping taxa.')
    tree1_lap = tree1.shear(taxa_lap)
    tree2_lap = tree2.shear(taxa_lap)
    return (tree1_lap, tree2_lap)


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
