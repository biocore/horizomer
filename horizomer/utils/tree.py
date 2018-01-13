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


def compare_topology(node1, node2):
    """Test whether the topologies of two trees are identical.

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
    Flipping child nodes does not change topology.
    Branch lengths are ignored.
    """
    if len(node1.children) != len(node2.children):
        return False
    elif len(node1.children) > 0:
        # check if sorted child names are identical
        children1 = sorted(node1.children, key=lambda x: x.name)
        children2 = sorted(node2.children, key=lambda x: x.name)
        taxa1 = [x.name for x in children1]
        taxa2 = [x.name for x in children2]
        if taxa1 != taxa2:
            return False

        # recursively test all children
        res = []
        for child1, child2 in zip(children1, children2):
            res.append(compare_topology(child1, child2))
        if res.count(True) < len(res):
            return False
    return True


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
        {TaxID : {'parent' : parent TaxID,
            'rank' : taxonomic rank,
            'name' : taxon name, empty if names_fp is None,
            'children' : set of (child TaxIDs)}}
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
        if pid not in taxdump:
            raise ValueError('TaxID %s is not defined.' % pid)
        if tid != pid:
            taxdump[pid]['children'].add(tid)

    return taxdump


def build_taxdump_tree(taxdump):
    """Build NCBI taxdump tree.

    Parameters
    ----------
    taxdump : dict of dict
        attributes of each TaxID, see read_taxdump

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
