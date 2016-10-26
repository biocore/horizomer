# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Reformat input files to format accepted by given HGT tool
=========================================================
"""

import click

from os import remove
from os.path import join
from skbio import TreeNode, TabularMSA, Sequence, Protein, DNA
from collections import OrderedDict


def join_trees(gene_tree,
               species_tree,
               output_tree_fp):
    """ Concatenate Newick trees into one file (species followed by gene).

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance for gene tree
    species_tree_fp: skbio.TreeNode
        TreeNode instance for species tree
    output_tree_fp: string
        file path to output species and gene tree

    See Also
    --------
    skbio.TreeNode
    """
    with open(output_tree_fp, 'w') as output_tree_f:
            output_tree_f.write(
                "%s\n%s\n" % (str(species_tree)[:-1], str(gene_tree)[:-1]))


def trim_gene_tree_leaves(gene_tree):
    """ Keep only string before first '_' delimiter in node ID.

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance

    See Also
    --------
    skbio.TreeNode

    Notes
    -----
    This function will keep only the word before the first '_' in the
    complete node ID. In ALF simulated sequences, the genes are labeled
    as "SPECIES_GENE". Most phylogenetic reconciliation tools
    require the associations between species leaves and gene leaves to
    be equal, therefore needing to remove the _GENENAME part in the gene
    tree.
    """
    for node in gene_tree.tips():
        node.name = node.name.split()[0]


def species_gene_mapping(gene_tree,
                         species_tree):
    """ Find the association between the leaves in species and gene trees.

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance for gene tree
    species_tree_fp: skbio.TreeNode
        TreeNode instance for species tree

    Returns
    -------
    mapping_leaves_t: OrderedDict
        Mapping between the species tree leaves and the gene tree leaves;
        species tips are the keys and gene tips are the values

    See Also
    --------
    skbio.TreeNode

    Notes
    -----
    Given the label format "SPECIES" for the species leaves and
    "SPECIES_GENE" in the gene leaves, report the associations between all
    species and gene leaves. Only one instance of the '_' delimiter is
    allowed in the gene leaves and this is used as a separator between the
    species name and the gene name.

    Ex.

    mapping = {"SE001":["SE001_1", "SE001_2"],
               "SE002":["SE002_1"]}
    """
    mapping_leaves = {}
    for node in species_tree.tips():
        if node.name not in mapping_leaves:
            mapping_leaves[node.name] = []
        else:
            raise ValueError(
                "Species tree leaves must be uniquely labeled: %s" % node.name)
    for node in gene_tree.tips():
        species, gene = node.name.split()
        if species in mapping_leaves:
            mapping_leaves[species].append("%s_%s" % (species, gene))
        else:
            raise ValueError(
                "Species %s does not exist in the species tree" % species)
    return OrderedDict(sorted(mapping_leaves.items(),
                       key=lambda x: x[1], reverse=True))


def remove_branch_lengths(tree):
    """ Set branch lengths to None.

    Parameters
    ----------
    tree: skbio.TreeNode
        TreeNode instance

    See Also
    --------
    skbio.TreeNode
    """
    for node in tree.postorder():
        node.length = None


def id_mapper(ids):
    mapping = {}
    for _id in ids:
        mapping[_id] = _id.split('/')[0]
    return mapping


def reformat_rangerdtl(gene_tree,
                       species_tree,
                       output_tree_fp):
    """ Reformat input trees to the format accepted by RANGER-DTL.

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance for gene tree
    species_tree_fp: skbio.TreeNode
        TreeNode instance for species tree
    output_tree_fp: string
        file path to output trees (species followed by gene)

    See Also
    --------
    skbio.TreeNode

    Notes
    -----
    The species name in the leaves of species and gene trees must be equal.
    For multiple genes from the same species, the format "SPECIES_GENE" is
    acceptable in the gene trees.

    """
    remove_branch_lengths(tree=gene_tree)
    remove_branch_lengths(tree=species_tree)
    join_trees(gene_tree,
               species_tree,
               output_tree_fp)


def reformat_trex(gene_tree,
                  species_tree,
                  output_tree_fp):
    """ Reformat input trees to the format accepted by T-REX.

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance for gene tree
    species_tree_fp: skbio.TreeNode
        TreeNode instance for species tree
    output_tree_fp: string
        file path to output trees (species followed by gene)

    See Also
    --------
    skbio.TreeNode

    Notes
    -----
    Binary trees only, leaves of species and gene trees must have equal
    names.
    """
    # trim gene tree leaves to exclude '_GENENAME' (if exists)
    trim_gene_tree_leaves(gene_tree)
    # join species and gene tree into one file
    join_trees(gene_tree,
               species_tree,
               output_tree_fp)


def reformat_riatahgt(gene_tree,
                      species_tree,
                      output_tree_fp):
    """ Reformat input trees to the format accepted by RIATA-HGT (PhyloNet).

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance for gene tree
    species_tree_fp: skbio.TreeNode
        TreeNode instance for species tree
    output_tree_fp: string
        file path to output trees (Nexus format)

    See Also
    --------
    skbio.TreeNode

    Notes
    -----
    Input to RIATA-HGT is a Nexus file. The number of leaves in the species
    and gene tree must be equal with the same naming.
    """
    nexus_file = """#NEXUS
BEGIN TREES;
Tree speciesTree = SPECIES_TREE
Tree geneTree = GENE_TREE
END;
BEGIN PHYLONET;
RIATAHGT speciesTree {geneTree};
END;
"""
    # trim gene tree leaves to exclude '_GENENAME' (if exists)
    trim_gene_tree_leaves(gene_tree)
    p = nexus_file.replace('SPECIES_TREE', str(species_tree)[:-1])
    p = p.replace('GENE_TREE', str(gene_tree)[:-1])
    with open(output_tree_fp, 'w') as output_tree_f:
        output_tree_f.write(p)


def reformat_jane4(gene_tree,
                   species_tree,
                   output_tree_fp):
    """ Reformat input trees to the format accepted by Jane4.

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance for gene tree
    species_tree_fp: skbio.TreeNode
        TreeNode instance for species tree
    output_tree_fp: string
        file path to output trees (Nexus format)

    See Also
    --------
    skbio.TreeNode

    Notes
    -----
    Input to Jane4 is a Nexus file, the trees cannot not contain branch
    lengths and the species/gene leaves mapping is required.
    """
    nexus_file = """#NEXUS
begin host;
tree host = SPECIES_TREE
endblock;
begin parasite;
tree parasite = GENE_TREE
endblock;
begin distribution;
Range MAPPING;
endblock;
"""
    # create a mapping between the species and gene tree leaves
    mapping_dict = species_gene_mapping(gene_tree=gene_tree,
                                        species_tree=species_tree)
    remove_branch_lengths(tree=gene_tree)
    remove_branch_lengths(tree=species_tree)
    mapping_str = ""
    for species in mapping_dict:
        for gene in mapping_dict[species]:
            mapping_str = "%s%s:%s, " % (mapping_str, gene, species)
    p = nexus_file.replace('SPECIES_TREE', str(species_tree))
    p = p.replace('GENE_TREE', str(gene_tree))
    p = p.replace('MAPPING', mapping_str[:-2])
    with open(output_tree_fp, 'w') as output_tree_f:
        output_tree_f.write(p)


def reformat_treepuzzle(gene_tree,
                        species_tree,
                        gene_msa_fa_fp,
                        output_tree_fp,
                        output_msa_phy_fp):
    """ Reformat input trees to the format accepted by Tree-Puzzle.

    Parameters
    ----------
    gene_tree: skbio.TreeNode
        TreeNode instance for gene tree
    species_tree_fp: skbio.TreeNode
        TreeNode instance for species tree
    gene_msa_fa_fp: string
        file path to gene alignments in FASTA format
    output_tree_fp: string
        file path to output trees (Nexus format)
    output_msa_phy_fp: string
        file path to output MSA in PHYLIP format

    See Also
    --------
    skbio.TreeNode
    """
    # remove the root branch length (output with ALF)
    for node in gene_tree.postorder():
        if node.is_root():
            node.length = None
    for node in species_tree.postorder():
        if node.is_root():
            node.length = None
    # trim gene tree leaves to exclude '_GENENAME' (if exists)
    trim_gene_tree_leaves(gene_tree)
    join_trees(gene_tree,
               species_tree,
               output_tree_fp)
    # trim FASTA sequence labels to exclude '/GENENAME' (if exists)
    msa_fa = TabularMSA.read(gene_msa_fa_fp, constructor=Protein)
    msa_fa.reassign_index(minter='id')
    mapping = id_mapper(msa_fa.index)
    msa_fa.reassign_index(mapping=mapping)
    msa_fa.write(output_msa_phy_fp, format='phylip')


def _merge_genbank_seqs(genbank_fp):
    """ Merge one to multiple sequences in a GenBank file into one.

    Parameters
    ----------
    genbank_fp: string
        file path to genome in GenBank format

    Returns
    -------
    tuple of (
        skbio.Sequence,
            Genome sequence, genes and metadata
        dict of { list of [ string, int, int, string ] }
            Gene name : translation, start, end, and strand
    )
    """
    loci = []
    nucl_seq = ''
    genes = {}
    nseq = 0  # number of nucleotide sequences
    with open(genbank_fp, 'r') as input_f:
        for line in input_f:
            if line.startswith('//'):
                nseq += 1
    abs_pos = 0  # absolute position in concantenated nucleotide sequence
    for i in range(nseq):
        gb = Sequence.read(genbank_fp, seq_num=i+1, format='genbank')
        locus_name = gb.metadata['LOCUS']['locus_name']
        size = gb.metadata['LOCUS']['size']
        loci.append([locus_name, size])
        nucl_seq += str(gb)
        for feature in gb.interval_metadata._intervals:
            m = feature.metadata
            if m['type'] == 'CDS' and 'protein_id' in m:
                protein_id = m['protein_id'].replace('\"', '')
                if protein_id not in genes:
                    translation = m['translation'].replace(' ', '') \
                        .replace('\"', '')
                    strand = m['strand']
                    start = feature.bounds[0][0] + abs_pos + 1
                    end = feature.bounds[0][1] + abs_pos
                    genes[protein_id] = [translation, start, end, strand]
        abs_pos += int(size)
    gb = DNA(nucl_seq)
    gb.metadata['LOCUS'] = {'locus_name': 'locus001', 'size': len(nucl_seq),
                            'unit': 'bp', 'shape': 'circular',
                            'division': 'CON', 'mol_type': 'DNA',
                            'date': '01-JAN-1900'}
    gb.metadata['id'] = 'locus001'
    gid = 1
    gb.interval_metadata._intervals = []
    for (gene, l) in sorted(genes.items(), key=lambda x: x[1][1]):
        location = str(l[1]) + '..' + str(l[2])
        if l[3] == '-':
            location = 'complement(' + location + ')'
        feature = {'type': 'gene', 'locus_tag': 'gene' + str(gid),
                   '__location': location}
        gb.interval_metadata.add([(l[1] - 1, l[2])], metadata=feature)
        feature = {'type': 'CDS', 'locus_tag': 'gene' + str(gid),
                   '__location': location, 'protein_id': gene,
                   'translation': l[0]}
        gb.interval_metadata.add([(l[1] - 1, l[2])], metadata=feature)
        gid += 1
    return (gb, genes)


def reformat_egid(genbank_fp,
                  output_dir):
    """ Reformat input genome to the formats accepted by EGID.

    Parameters
    ----------
    genbank_fp: string
        file path to genome in GenBank format
    output_dir: string
        output directory path

    Notes
    -----
    Input to EGID are five obsolete NCBI standard files: gbk, fna, faa, ffn
    and ptt.
    """
    (gb, genes) = _merge_genbank_seqs(genbank_fp)
    DNA.write(gb, join(output_dir, 'id.fna'), format='fasta')
    DNA.write(gb, join(output_dir, 'id.gbk'), format='genbank')
    nucl_seq = str(gb)
    output_f = {}
    for x in ('faa', 'ffn', 'ptt'):
        output_f[x] = open(join(output_dir, 'id.' + x), 'w')
    output_f['ptt'].write('locus001\n' + str(len(genes)) + ' proteins\n')
    fields = ('Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym',
              'Code', 'COG', 'Product')
    output_f['ptt'].write('\t'.join(fields) + '\n')
    gid = 1
    for (gene, l) in sorted(genes.items(), key=lambda x: x[1][1]):
        output_f['faa'].write('>' + gene + '\n' + l[0] + '\n')
        output_f['ptt'].write(str(l[1]) + '..' + str(l[2]) + '\t' +
                              l[3] + '\t' + str(len(l[0])) + '\t' +
                              str(gid) + '\t-\tgene' + str(gid) +
                              '\t-\t-\t-\n')
        if l[3] == '+':
            output_f['ffn'].write('>locus001:' + str(l[1]) + '-' +
                                  str(l[2]) + '\n' +
                                  nucl_seq[l[1]-1:l[2]] + '\n')
        else:
            rc_seq = str(DNA(nucl_seq[l[1]-1:l[2]]).reverse_complement())
            output_f['ffn'].write('>locus001:c' + str(l[2]) + '-' +
                                  str(l[1]) + '\n' + rc_seq + '\n')
        gid += 1
    for x in output_f:
        output_f[x].close()


def reformat_genemark(genbank_fp,
                      output_dir):
    """ Reformat input genome to the formats accepted by GeneMark.

    Parameters
    ----------
    genbank_fp: string
        file path to genome in GenBank format
    output_dir: string
        output directory path

    Notes
    -----
    GeneMark's acceptable input file format is FASTA (genome sequence).
    """
    gb = _merge_genbank_seqs(genbank_fp)[0]
    DNA.write(gb, join(output_dir, 'id.fna'), format='fasta')
    DNA.write(gb, join(output_dir, 'id.gbk'), format='genbank')


@click.command()
@click.option('--gene-tree-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Gene tree in Newick format')
@click.option('--species-tree-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Species tree in Newick format')
@click.option('--gene-msa-fa-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='MSA of genes in FASTA format')
@click.option('--genbank-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Genome in GenBank format')
@click.option('--output-tree-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output formatted species and gene tree')
@click.option('--output-msa-phy-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output MSA in PHYLIP format')
@click.option('--output-dir', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output directory path')
@click.option('--method', required=True,
              type=click.Choice(['trex', 'ranger-dtl',
                                 'riata-hgt', 'consel',
                                 'darkhorse', 'hgtector',
                                 'genemark', 'egid',
                                 'distance-method', 'jane4',
                                 'wn-svm', 'tree-puzzle']),
              help='The method to be used for HGT detection')
def _main(gene_tree_fp,
          species_tree_fp,
          gene_msa_fa_fp,
          genbank_fp,
          output_tree_fp,
          output_msa_phy_fp,
          output_dir,
          method):
    """ Reformat input files to format accepted by various HGT tools.

    For phylogenetic methods, a species tree and a gene tree are mandatory.
    Species tree can be multifurcating, however will be converted to
    bifurcating trees for software that require them. Leaf labels of
    species tree and gene tree must match, however the label
    SPECIES_GENE is acceptable for multiple genes in the gene
    tree. Leaf labels must also be at most 10 characters long (for
    PHYLIP manipulations).

    For compositional methods, a GenBank file containing both the genome
    sequence and the coordinates of its gene regions is required. Draft
    genomes (multiple sequences) are acceptable.
    """

    # add function to check where tree is multifurcating and the labeling
    # is correct
    gene_tree = TreeNode.read(gene_tree_fp, format='newick') \
        if gene_tree_fp is not None else None
    species_tree = TreeNode.read(species_tree_fp, format='newick') \
        if species_tree_fp is not None else None

    if method == 'ranger-dtl':
        reformat_rangerdtl(
            gene_tree=gene_tree,
            species_tree=species_tree,
            output_tree_fp=output_tree_fp)
    elif method == 'trex':
        reformat_trex(
            gene_tree=gene_tree,
            species_tree=species_tree,
            output_tree_fp=output_tree_fp)
    elif method == 'riata-hgt':
        reformat_riatahgt(
            gene_tree=gene_tree,
            species_tree=species_tree,
            output_tree_fp=output_tree_fp)
    elif method == 'jane4':
        reformat_jane4(
            gene_tree=gene_tree,
            species_tree=species_tree,
            output_tree_fp=output_tree_fp)
    elif method == 'tree-puzzle':
        reformat_treepuzzle(
            gene_tree=gene_tree,
            species_tree=species_tree,
            gene_msa_fa_fp=gene_msa_fa_fp,
            output_tree_fp=output_tree_fp,
            output_msa_phy_fp=output_msa_phy_fp)
    elif method == 'egid':
        reformat_egid(
            genbank_fp=genbank_fp,
            output_dir=output_dir)
    elif method == 'genemark':
        reformat_genemark(
            genbank_fp=genbank_fp,
            output_dir=output_dir)


if __name__ == "__main__":
    _main()
