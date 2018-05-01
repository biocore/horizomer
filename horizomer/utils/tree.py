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


def order_nodes(tree, increase=True):
    """Rotate internal nodes of a tree so that child nodes are ordered by the
    number of descendants.

    Parameters
    ----------
    tree : skbio.TreeNode
        tree to order
    increase : bool, optional
        order nodes in increasing (True) or decreasing (False) order

    Returns
    -------
    skbio.TreeNode
        resulting ordered tree

    See Also
    --------
    is_ordered

    Examples
    --------
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(['(((a,b),(c,d,e)),((f,g),h));'])
    >>> print(tree)
    (((a,b),(c,d,e)),((f,g),h));
    <BLANKLINE>
    >>> tree_ordered = order_nodes(tree, False)
    >>> print(tree_ordered)
    ((h,(f,g)),((a,b),(c,d,e)));
    <BLANKLINE>
    """
    res = tree.copy()
    for node in res.postorder():
        if node.is_tip():
            node.n = 1
        else:
            node.n = sum(x.n for x in node.children)
    for node in res.postorder():
        if not node.is_tip():
            children = node.children
            node.children = []
            for child in sorted(children, key=lambda x: x.n, reverse=increase):
                node.append(child)
    for node in res.postorder():
        delattr(node, 'n')
    return res


def is_ordered(tree, increase=True):
    """Returns `True` if the tree is ordered.

    Parameters
    ----------
    tree : skbio.TreeNode
        tree to check ordering
    increase : bool, optional
        check if nodes in increasing (True) or decreasing (False) order

    Returns
    -------
    bool
        `True` if the tree is ordered

    See Also
    --------
    order_nodes

    Examples
    --------
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(['((a,b)c,d)e;'])
    >>> is_ordered(tree)
    True
    """
    tcopy = tree.copy()
    for node in tcopy.postorder():
        if node.is_tip():
            node.n = 1
        else:
            node.n = sum(x.n for x in node.children)

    p = tcopy.root()
    prev = p.n
    for node in tcopy.levelorder():
        s = [x for x in p.siblings()]
        if node in s:
            cur = node.n
            if prev < cur if increase else prev > cur:
                return False
            prev = cur
        else:
            p = node
            prev = p.n
    return True


def cladistic(tree, taxa):
    """Determines the cladistic property of the given taxon set.

    Parameters
    ----------
    tree : skbio.TreeNode
        tree for taxa comparison
    taxa : iterable of str
        taxon names

    Returns
    -------
    str
        'uni' if input taxon is a single tip in given tree
        'mono' if input taxa are monophyletic in given tree
        'poly' if input taxa are polyphyletic in given tree

    Notes
    -----
    In the following tree example:
                                  /-a
                        /--------|
                       |          \-b
              /--------|
             |         |          /-c
             |         |         |
             |          \--------|--d
    ---------|                   |
             |                    \-e
             |
             |                    /-f
             |          /--------|
              \--------|          \-g
                       |
                        \-h
    ['a'] returns 'uni'
    ['c', 'd', 'e'] returns 'mono'
    ['a', 'c', 'f'] returns 'poly'
    ['f', 'h'] returns 'poly'
    Paraphyly, which is programmably indistinguishable from polyphyly, returns
    poly here.

    Raises
    ------
    ValueError
        if one or more taxon names are not present in the tree

    Examples
    --------
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(['((a,b)c,d)e;'])
    >>> cladistic(tree, ['a'])
    'uni'
    >>> cladistic(tree, ['a', 'b'])
    'mono'
    >>> cladistic(tree, ['a', 'd'])
    'poly'
    """
    tips = []
    taxa = set(taxa)
    for tip in tree.tips():
        if tip.name in taxa:
            tips.append(tip)
    n = len(taxa)
    if len(tips) < n:
        raise ValueError('Taxa not found in the tree.')
    return ('uni' if n == 1 else
            ('mono' if len(tree.lca(tips).subset()) == n else
             'poly'))


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

    Examples
    --------
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(['((a,b)99,(c,d):1.0);'])
    >>> support(tree.lca(['a', 'b']))
    99.0
    >>> support(tree.lca(['c', 'd'])) is None
    True
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

    Raises
    ------
    ValueError
        if input node is root

    Examples
    --------
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(['((c:2.0,d:3.0)a:1.0,(e:2.0,f:1.0)b:2.0);'])
    >>> unpack(tree.find('b'))
    >>> print(tree)
    ((c:2.0,d:3.0)a:1.0,e:4.0,f:3.0);
    <BLANKLINE>
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

    Raises
    ------
    ValueError
        if taxon is empty
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


def unpack_by_func(tree, func):
    """Unpack internal nodes that meet certain criteria.

    Parameters
    ----------
    tree : skbio.TreeNode
        tree to search for nodes to unpack
    func : function
        a function that accepts a TreeNode and returns `True` or `False`,
        where `True` indicates the node is to be unpacked

    Returns
    -------
    skbio.TreeNode
        resulting tree with nodes meeting criteria unpacked

    Examples
    --------
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(['((c:2,d:3)a:1,(e:1,f:2)b:2);'])
    >>> tree_unpacked = unpack_by_func(tree, lambda x: x.length <= 1)
    >>> print(tree_unpacked)
    ((e:1.0,f:2.0)b:2.0,c:3.0,d:4.0);
    <BLANKLINE>
    >>> tree = TreeNode.read(['(((a,b)85,(c,d)78)75,(e,(f,g)64)80);'])
    >>> tree_unpacked = unpack_by_func(tree, lambda x: support(x) < 75)
    >>> print(tree_unpacked)
    (((a,b)85,(c,d)78)75,(e,f,g)80);
    <BLANKLINE>
    """
    tcopy = tree.copy()
    nodes_to_unpack = []
    for node in tcopy.non_tips():
        if func(node):
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
