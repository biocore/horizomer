# this script generates a tree of a set genomes based on their NCBI taxonomy
# (taxdump) annotations
# usage: python me.py nodes.dmp genome2taxid.tsv output.nwk

import sys
from skbio import TreeNode

# read NCBI taxonomy hierarchies (nodes.dmp)
taxdump = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        # format of nodes.dmp:
        # taxid <tab> | <tab> parent taxid <tab> | <tab> more info...
        x = line.rstrip('\r\n').split('\t|\t')
        taxdump[x[0]] = {'parent': x[1], 'children': set()}

# identify root of the tree and children of all nodes
root = None
for tid in taxdump:
    pid = taxdump[tid]['parent']
    # at root, the taxid and its parental taxid are both 1
    if tid == pid:
        root = tid
    elif pid in taxdump:
        taxdump[pid]['children'].add(tid)
    else:
        raise ValueError('TaxID %s is not defined.' % pid)

# create the tree from root
tree = TreeNode(root)


# iteratively attach child nodes to parent node
def iter_node(node):
    for cid in taxdump[node.name]['children']:
        child = TreeNode(cid)
        node.extend([child])
        iter_node(child)


iter_node(tree)

print('The taxdump tree has %d nodes, including %d tips.'
      % (tree.count(), tree.count(tips=True)))


# read genome ID to TaxID translation table
# format: genome ID <tab> TaxID
# note: multiple genome IDs may correspond to one TaxID
with open(sys.argv[2], 'r') as f:
    g2tid = dict(x.split('\t') for x in f.read().splitlines())
tid2gs = {}  # TaxID to genome ID dictionary
for g, tid in g2tid.items():
    if tid not in tid2gs:
        tid2gs[tid] = set([g])
    else:
        tid2gs[tid].add(g)
tids = set(tid2gs)  # TaxIDs to be retained
print('The dataset has %d genomes assigned to %d TaxIDs.'
      % (len(g2tid), len(tids)))


# shrink the taxonomy tree to only contain tips present in the dataset
# logic: recursively remove tips which are not in the dataset, until no
# more tips can be removed.
while True:
    tips_to_remove = set()
    for tip in tree.tips():
        if tip.name not in tids:
            tips_to_remove.add(tip)
    n = len(tips_to_remove)
    if n == 0:
        break
    tree.remove_deleted(lambda x: x in tips_to_remove)

print('The taxdump tree was shrinked to %d nodes, including %d tips, all of '
      'which are in the dataset.' % (tree.count(), tree.count(tips=True)))

# collapse internal nodes with single child and not in the dataset
# modified from scikit-bio's `prune` function
nodes_to_remove = []
for node in tree.traverse():
    if len(node.children) == 1 and node.name not in tids:
        nodes_to_remove.append(node)
for node in nodes_to_remove:
    if node.is_root():
        tree = node.children[0]
    else:
        node.parent.append(node.children[0])
        node.parent.remove(node)
print('%d single-child internal nodes were removed. The tree now has %d nodes.'
      % (len(nodes_to_remove), tree.count()))

# for internal nodes that are in the dataset, raise their childrens to be
# siblings of them, so that themselves become tips
nodes_to_collapse = []
for node in tree.non_tips():
    if node.name in tids:
        nodes_to_collapse.append(node)
for node in nodes_to_collapse:
    node.parent.extend(node.children)
    for child in node.children:
        node.remove(child)
print('For %d internal nodes that are in the dataset, their child nodes were '
      'raised to be their siblings.' % len(nodes_to_collapse))

# validate that the tree transformation is complete
for tip in tree.tips():
    if tip.name not in tids:
        raise ValueError('Error: Failed to remove tip %s.' % tip.name)
for node in tree.non_tips():
    if node.name in tids:
        raise ValueError('Error: Failed to convert %s from an internal node '
                         'to a tip.' % node.name)
n = tree.count(tips=True) - len(tids)
if n > 0:
    raise ValueError('Error: Failed to retain %d TaxIDs in the dataset as '
                     'tips of the tree.' % n)
print('The tree now has %d nodes, including %d tips. The tips and the TaxIDs '
      'in the dataset match.' % (tree.count(), tree.count(tips=True)))

# translate TaxIDs into genome IDs
# if multiple genome IDs correspond to one TaxID, make them a polytomy
tips = [x for x in tree.tips()]
for tip in tips:
    tip.parent.extend([TreeNode(x) for x in tid2gs[tip.name]])
    tip.parent.remove(tip)
print('The TaxIDs at the tips were translated into genome IDs.')

# validate that the tree translation is complete
for tip in tree.tips():
    if tip.name not in g2tid:
        raise ValueError('Error: Failed to translate tip %s.' % tip.name)
n = tree.count(tips=True) - len(g2tid)
if n > 0:
    raise ValueError('Error: Failed to retain %d genome IDs as tips of the '
                     'tree.' % n)
print('The tree now has %d nodes, including %d tips. The tips and the genome '
      'IDs in the dataset match.' % (tree.count(), tree.count(tips=True)))

# output tree in Newick format
tree.write(sys.argv[3])
