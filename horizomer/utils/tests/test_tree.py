#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The Horizomer Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join, dirname, realpath
from skbio import TreeNode

from horizomer.utils.tree import (
    support, collapse, has_duplicates, compare_topology, intersect_trees,
    collapse_short_branches, collapse_low_support_nodes,
    read_taxdump, build_taxdump_tree)


class TreeTests(TestCase):

    def setUp(self):
        """ Set up working directory and test files
        """
        # test output can be written to this directory
        self.working_dir = mkdtemp()

        # test data directory
        datadir = join(dirname(realpath(__file__)), 'data')

        # test data files
        self.nodes_fp = join(datadir, 'nodes.dmp')
        self.names_fp = join(datadir, 'names.dmp')

    def tearDown(self):
        # there isn't any file to remove at the moment
        # but in the future there will be
        rmtree(self.working_dir)

    def test_support(self):
        """Test getting support value of a node."""
        # test nodes with support alone as label
        tree = TreeNode.read(['((a,b)75,(c,d)90);'])
        node1, node2 = tree.children
        self.assertEqual(support(node1), 75.0)
        self.assertEqual(support(node2), 90.0)

        # test nodes with support and branch length
        tree = TreeNode.read(['((a,b)0.85:1.23,(c,d)0.95:4.56);'])
        node1, node2 = tree.children
        self.assertEqual(support(node1), 0.85)
        self.assertEqual(support(node2), 0.95)

        # test nodes with support and extra label (not a common scenario but
        # can happen)
        tree = TreeNode.read(['((a,b)\'80:X\',(c,d)\'60:Y\');'])
        node1, node2 = tree.children
        self.assertEqual(support(node1), 80.0)
        self.assertEqual(support(node2), 60.0)

        # test nodes without label, with non-numeric label, and with branch
        # length only
        tree = TreeNode.read(['((a,b),(c,d)x,(e,f):1.0);'])
        for node in tree.children:
            self.assertIsNone(support(node))

    def test_collapse(self):
        """Test collapsing an internal node."""
        # test collapsing a node without branch length
        tree = TreeNode.read(['((c,d)a,(e,f)b);'])
        collapse(tree.find('b'))
        exp = '((c,d)a,e,f);\n'
        self.assertEqual(str(tree), exp)

        # test collapsing a node with branch length
        tree = TreeNode.read(['((c:2.0,d:3.0)a:1.0,(e:2.0,f:1.0)b:2.0);'])
        collapse(tree.find('b'))
        exp = '((c:2.0,d:3.0)a:1.0,e:4.0,f:3.0);'
        self.assertEqual(str(tree).rstrip(), exp)

        # test attempting to collapse root
        tree = TreeNode.read(['((d,e)b,(f,g)c)a;'])
        msg = 'Cannot collapse root.'
        with self.assertRaisesRegex(ValueError, msg):
            collapse(tree.find('a'))

    def test_has_duplicates(self):
        """Test checking for duplicated taxa."""
        # test tree without duplicates
        tree = '((a,b),(c,d));'
        obs = has_duplicates(TreeNode.read([tree]))
        self.assertFalse(obs)

        # test tree with duplicates
        tree = '((a,a),(c,a));'
        obs = has_duplicates(TreeNode.read([tree]))
        self.assertTrue(obs)

        tree = '((1,(2,x)),4,(5,(6,x,8)));'
        obs = has_duplicates(TreeNode.read([tree]))
        self.assertTrue(obs)

        # test tree with empty taxon names (not a common scenario but can
        # happen)
        tree = '((1,(2,,)),4,(5,(6,,8)));'
        msg = 'Empty taxon name\(s\) found.'
        with self.assertRaisesRegex(ValueError, msg):
            has_duplicates(TreeNode.read([tree]))

    def test_compare_topology(self):
        """Test comparing topologies of two trees."""
        # test identical Newick strings
        tree1 = '(a,b)c;'
        tree2 = '(a,b)c;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertTrue(obs)

        # test identical topologies with different branch lengths
        tree1 = '(a:1,b:2)c:3;'
        tree2 = '(a:3,b:2)c:1;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertTrue(obs)

        # test identical topologies with flipped child nodes
        tree1 = '(a,b)c;'
        tree2 = '(b,a)c;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertTrue(obs)

        tree1 = '((4,5)2,(6,7,8)3)1;'
        tree2 = '((8,7,6)3,(5,4)2)1;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertTrue(obs)

        tree1 = '(((9,10)4,(11,12,13)5)2,((14)6,(15,16,17,18)7,(19,20)8)3)1;'
        tree2 = '(((15,16,17,18)7,(14)6,(20,19)8)3,((12,13,11)5,(10,9)4)2)1;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertTrue(obs)

        # test different topologies
        tree1 = '(a,b)c;'
        tree2 = '(a,c)b;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertFalse(obs)

        tree1 = '((4,5)2,(6,7,8)3)1;'
        tree2 = '((4,5)3,(6,7,8)2)1;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertFalse(obs)

        tree1 = '((4,5)2,(6,7,8)3)1;'
        tree2 = '(((4,1)8)7,(6,3)2)5;'
        obs = compare_topology(TreeNode.read([tree1]), TreeNode.read([tree2]))
        self.assertFalse(obs)

    def test_intersect_trees(self):
        """Test intersecting two trees."""
        # test trees with identical taxa
        tree1 = '((a,b),(c,d));'
        tree2 = '(a,(b,c,d));'
        obs = intersect_trees(TreeNode.read([tree1]), TreeNode.read([tree2]))
        exp = (TreeNode.read([tree1]), TreeNode.read([tree2]))
        for i in range(2):
            self.assertEqual(obs[i].compare_subsets(exp[i]), 0.0)

        # test trees with partially different taxa
        tree1 = '((a,b),(c,d));'
        tree2 = '((a,b),(c,e));'
        obs = intersect_trees(TreeNode.read([tree1]), TreeNode.read([tree2]))
        tree1_lap = '((a,b),c);'
        tree2_lap = '((a,b),e);'
        exp = (TreeNode.read([tree1_lap]), TreeNode.read([tree2_lap]))
        for i in range(2):
            self.assertEqual(obs[i].compare_subsets(exp[i]), 0.0)

        tree1 = '(((a,b),(c,d)),((e,f,g),h));'
        tree2 = '(a,((b,x),(d,y,(f,g,h))));'
        obs = intersect_trees(TreeNode.read([tree1]), TreeNode.read([tree2]))
        tree1_lap = '(((a,b),d),((f,g),h));'
        tree2_lap = '(a,(b,(d,(f,g,h))));'
        exp = (TreeNode.read([tree1_lap]), TreeNode.read([tree2_lap]))
        for i in range(2):
            self.assertEqual(obs[i].compare_subsets(exp[i]), 0.0)

        # test trees with completely different taxa
        tree1 = '((a,b),(c,d));'
        tree2 = '((e,f),(g,h));'
        msg = 'Trees have no overlapping taxa.'
        with self.assertRaisesRegex(ValueError, msg):
            intersect_trees(TreeNode.read([tree1]), TreeNode.read([tree2]))

        # test trees with duplicated taxa
        tree1 = '((a,b),(c,d));'
        tree2 = '((a,a),(b,c));'
        msg = 'Either tree has duplicated taxa.'
        with self.assertRaisesRegex(ValueError, msg):
            intersect_trees(TreeNode.read([tree1]), TreeNode.read([tree2]))

    def test_collapse_short_branches(self):
        """Test collapsing short branches."""
        # collapse internal nodes with branch length < 1.01
        # (will collapse node 'a', but not tip 'e')
        # (will add the branch length of 'a' to its child nodes 'c' and 'd')
        tree = TreeNode.read(['((c:2,d:3)a:1,(e:1,f:2)b:2);'])
        obs = str(collapse_short_branches(tree, 1.01)).rstrip()
        exp = '((e:1.0,f:2.0)b:2.0,c:3.0,d:4.0);'
        self.assertEqual(obs, exp)

        # collapse internal nodes with branch length < 2.01
        # (will collapse both 'a' and 'b')
        obs = str(collapse_short_branches(tree, 2.01)).rstrip()
        exp = '(c:3.0,d:4.0,e:3.0,f:4.0);'
        self.assertEqual(obs, exp)

        # collapse two nested nodes 'a' and 'c' simultaneously
        tree = TreeNode.read(['(((e:3,f:2)c:1,d:3)a:1,b:4);'])
        obs = str(collapse_short_branches(tree, 2.01)).rstrip()
        exp = '(b:4.0,d:4.0,e:5.0,f:4.0);'
        self.assertEqual(obs, exp)

        # collapse internal node 'a' which has no branch length
        tree = TreeNode.read(['((c:2,d:3)a,(e:1,f:2)b:2);'])
        obs = str(collapse_short_branches(tree, 1.01)).rstrip()
        exp = '((e:1.0,f:2.0)b:2.0,c:2.0,d:3.0);'
        self.assertEqual(obs, exp)

        # test a complicated scenario (collapsing nodes 'g', 'h' and 'm')
        tree = TreeNode.read(['(((a:1.04,b:2.32,c:1.44)d:3.20,'
                              '(e:3.91,f:2.47)g:1.21)h:1.75,'
                              '(i:4.14,(j:2.06,k:1.58)l:3.32)m:0.77);'])
        obs = str(collapse_short_branches(tree, 2.01)).rstrip()
        exp = ('((a:1.04,b:2.32,c:1.44)d:4.95,e:6.87,f:5.43,i:4.91,'
               '(j:2.06,k:1.58)l:4.09);')
        self.assertEqual(obs, exp)

    def test_collapse_low_support_nodes(self):
        """Test collapsing low-support nodes."""
        # collapse nodes with support < 75
        tree = TreeNode.read(['(((a,b)85,(c,d)78)75,(e,(f,g)64)80);'])
        obs = str(collapse_low_support_nodes(tree, 75)).rstrip()
        exp = '(((a,b)85,(c,d)78)75,(e,f,g)80);'
        self.assertEqual(obs, exp)

        # collapse nodes with support < 85
        obs = str(collapse_low_support_nodes(tree, 85)).rstrip()
        exp = '((a,b)85,c,d,e,f,g);'
        self.assertEqual(obs, exp)

        # collapse nodes with support < 0.95
        tree = TreeNode.read(['(((a,b)0.97,(c,d)0.98)1.0,(e,(f,g)0.88)0.96);'])
        obs = str(collapse_low_support_nodes(tree, 0.95)).rstrip()
        exp = '(((a,b)0.97,(c,d)0.98)1.0,(e,f,g)0.96);'
        self.assertEqual(obs, exp)

        # test a case where there are branch lengths, none support values add
        # node labels
        tree = TreeNode.read(['(((a:1.02,b:0.33)85:0.12,(c:0.86,d:2.23)'
                              '70:3.02)75:0.95,(e:1.43,(f:1.69,g:1.92)64:0.20)'
                              'node:0.35)root;'])
        obs = str(collapse_low_support_nodes(tree, 75)).rstrip()
        exp = ('(((a:1.02,b:0.33)85:0.12,c:3.88,d:5.25)75:0.95,'
               '(e:1.43,f:1.89,g:2.12)node:0.35)root;')
        self.assertEqual(obs, exp)

    def test_read_taxdump(self):
        """Test reading NCBI taxdump."""
        obs = read_taxdump(self.nodes_fp)
        exp = {
            '1': {'parent': '1', 'rank': 'order',
                  'children': set(['2', '3'])},
            '2': {'parent': '1', 'rank': 'family',
                  'children': set(['4', '5'])},
            '3': {'parent': '1', 'rank': 'family',
                  'children': set(['6', '7', '8'])},
            '4': {'parent': '2', 'rank': 'genus',
                  'children': set(['9', '10'])},
            '5': {'parent': '2', 'rank': 'genus',
                  'children': set(['11', '12', '13'])},
            '6': {'parent': '3', 'rank': 'genus',
                  'children': set(['14'])},
            '7': {'parent': '3', 'rank': 'genus',
                  'children': set(['15', '16', '17', '18'])},
            '8': {'parent': '3', 'rank': 'genus',
                  'children': set(['19', '20'])},
            '9': {'parent': '4', 'rank': 'species', 'children': set()},
            '10': {'parent': '4', 'rank': 'species', 'children': set()},
            '11': {'parent': '5', 'rank': 'species', 'children': set()},
            '12': {'parent': '5', 'rank': 'species', 'children': set()},
            '13': {'parent': '5', 'rank': 'species', 'children': set()},
            '14': {'parent': '6', 'rank': 'species', 'children': set()},
            '15': {'parent': '7', 'rank': 'species', 'children': set()},
            '16': {'parent': '7', 'rank': 'species', 'children': set()},
            '17': {'parent': '7', 'rank': 'species', 'children': set()},
            '18': {'parent': '7', 'rank': 'species', 'children': set()},
            '19': {'parent': '8', 'rank': 'species', 'children': set()},
            '20': {'parent': '8', 'rank': 'species', 'children': set()}
        }
        for tid in exp:
            exp[tid]['name'] = ''
        self.assertDictEqual(obs, exp)

        obs = read_taxdump(self.nodes_fp, self.names_fp)
        name_dict = {
            '1': 'root', '2': 'Eukaryota', '3': 'Bacteria', '4': 'Plantae',
            '5': 'Animalia', '6': 'Bacteroidetes', '7': 'Proteobacteria',
            '8': 'Firmicutes', '9': 'Gymnosperms', '10': 'Angiosperms',
            '11': 'Chordata', '12': 'Arthropoda', '13': 'Mollusca',
            '14': 'Prevotella', '15': 'Escherichia', '16': 'Vibrio',
            '17': 'Rhizobium', '18': 'Helicobacter', '19': 'Bacillus',
            '20': 'Clostridia'
        }
        for tid in name_dict:
            exp[tid]['name'] = name_dict[tid]
        self.assertDictEqual(obs, exp)

    def test_build_taxdump_tree(self):
        """Test building NCBI taxdump tree."""
        taxdump = read_taxdump(self.nodes_fp)
        obs = build_taxdump_tree(taxdump)
        exp = TreeNode.read(['(((9,10)4,(11,12,13)5)2,((14)6,(15,16,17,18)7,'
                             '(19,20)8)3)1;'])
        self.assertTrue(compare_topology(obs, exp))


if __name__ == '__main__':
    main()
