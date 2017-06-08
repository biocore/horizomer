## Installation

(updated on 2/9/2017)

The Horizomer pipeline relies on multiple external applications. We are working to simplify the installation process. Please follow the instructions in this document.

### Operating system

A Unix-like operating system (e.g., Linux or Mac OS) is required to run Horizomer.

Note: The pipeline has been tested on CentOS 6.6, Ubuntu 16.04, and Mac OS 10.11.

### Recommended installation procedures

Have [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://anaconda.org/) installed in the system. The Python 3 version is preferable.

Create a conda environment and install required libraries:
```
conda create -n horizomer python=3.5
source activate horizomer
conda install click biopython
conda install -c biocore scikit-bio
```

Install third-party applications using conda:
```
conda install -c bioconda blast diamond fasttree mafft mcl muscle prodigal=2.6.2 raxml t_coffee trimal
```

Exit the conda environment when done.
```
source deactivate
```

For applications that require Python 2 (specifically: **PhyloPhlAn** and **OrthoFinder**), create another conda environment and install them:
```
conda create -n horizomer_py2 -c bioconda python=2.7 biopython orthofinder
```

The remaining applications have to be installed manually. Please refer to the table below for details.

### List of required external applications

(Note: those with "how to install" = "conda" are already installed if you have run the aforementioned commands.)

| Name | Tested Version | Purpose | How to install | PMID |
| --- | --- | --- | --- | --- |
| [T-REX](http://www.trex.uqam.ca/index.php?action=hgt&project=trex) | [3.22](http://www.trex.uqam.ca/download/hgt-detection_3.22.zip) | HGT detection | download & compile | [20525630](https://www.ncbi.nlm.nih.gov/pubmed/20525630) |
| [RANGER-DTL](http://compbio.engr.uconn.edu/software/RANGER-DTL/) | [2.0](http://compbio.engr.uconn.edu/software/RANGER-DTL/Linux.zip) | HGT detection | download binary | [22689773](https://www.ncbi.nlm.nih.gov/pubmed/22689773) |
| [PhyloNet](https://bioinfocs.rice.edu/phylonet) | [3.6.1](https://bioinfocs.rice.edu/sites/g/files/bxs266/f/kcfinder/files/PhyloNet_3.6.1.jar) | HGT detection | download binary | [18662388](https://www.ncbi.nlm.nih.gov/pubmed/18662388) |
| [Jane](https://www.cs.hmc.edu/~hadas/jane/index.html) | [4.01](https://www.cs.hmc.edu/~hadas/jane/form.html) | HGT detection | download binary (!license!) | [20181081](https://www.ncbi.nlm.nih.gov/pubmed/20181081) |
| [TREE-PUZZLE](http://www.tree-puzzle.de/) | [5.3.rc16](http://www.tree-puzzle.de/tree-puzzle-5.3.rc16-linux.tar.gz) | HGT detection | download & compile | [11934758](https://www.ncbi.nlm.nih.gov/pubmed/11934758) |
| [CONSEL](http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/) | [0.20](http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/pub/cnsls020.tgz) | HGT detection | download | [11751242](https://www.ncbi.nlm.nih.gov/pubmed/11751242) |
| [DarkHorse](http://darkhorse.ucsd.edu/) | [1.5 rev170](http://darkhorse.ucsd.edu/DarkHorse-1.5_rev170.tar.gz) | HGT detection | download & install | [17274820](https://www.ncbi.nlm.nih.gov/pubmed/17274820) |
| [HGTector](https://github.com/DittmarLab/HGTector) | [0.2.1](https://github.com/DittmarLab/HGTector/archive/wgshgt.zip) | HGT detection | git clone | [25159222](https://www.ncbi.nlm.nih.gov/pubmed/25159222) |
| [EGID](http://www5.esu.edu/cpsc/bioinfo/software/EGID/) | [1.0](http://www5.esu.edu/cpsc/bioinfo/software/EGID/EGID_1.0.tar.gz) | HGT detection | download | [22355228](https://www.ncbi.nlm.nih.gov/pubmed/22355228) |
| [GeneMarkS](http://exon.gatech.edu/GeneMark/) | [4.30](http://exon.gatech.edu/GeneMark/license_download.cgi) | HGT detection | download binary (!license!) | [9461475](https://www.ncbi.nlm.nih.gov/pubmed/9461475) |
||
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) | [1.1.4](https://github.com/davidemms/OrthoFinder/releases/download/1.1.4/OrthoFinder-1.1.4.tar.gz) | orthology identification | conda | [26243257](https://www.ncbi.nlm.nih.gov/pubmed/26243257) |
| [Phylomizer](https://github.com/Gabaldonlab/phylomizer) | [9/12/2016](https://github.com/Gabaldonlab/phylomizer/commit/e427a04b3d62bbac4d760fef975f6bdf5aeed44a) | gene family tree building | git clone | NA |
| [PhyloPhlAn](https://github.com/davidemms/OrthoFinder) | [1/25/2017](https://bitbucket.org/nsegata/phylophlan/commits/2c0e61ad820b8ff732837e98f7843afbb7ec1cda) | genome tree building | hg clone | [23942190](https://www.ncbi.nlm.nih.gov/pubmed/23942190) |
||
| [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) | [2.5.0](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/ncbi-blast-2.5.0+-src.tar.gz) | sequence similarity search | conda | [2231712](https://www.ncbi.nlm.nih.gov/pubmed/2231712) |
| [DIAMOND](https://ab.inf.uni-tuebingen.de/software/diamond) | [0.8.28](https://github.com/bbuchfink/diamond/archive/v0.8.28.tar.gz) | sequence similarity search | conda | [25402007](https://www.ncbi.nlm.nih.gov/pubmed/25402007) |
| [USEARCH](http://www.drive5.com/usearch/) | [9.1.13](http://www.drive5.com/usearch/download.html) | sequence similarity search | download binary (!license!) | [20709691](https://www.ncbi.nlm.nih.gov/pubmed/20709691) |
||
| [MAFFT](http://mafft.cbrc.jp/alignment/software/) | [7.305](http://mafft.cbrc.jp/alignment/software/mafft-7.305-with-extensions-src.tgz)  | sequence alignment | conda | [12136088](https://www.ncbi.nlm.nih.gov/pubmed/12136088) |
| [MUSCLE](http://drive5.com/muscle/) | [3.8.31](http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz) | sequence alignment | conda | [15034147](https://www.ncbi.nlm.nih.gov/pubmed/15034147) |
| [KAlign](http://www.ebi.ac.uk/Tools/msa/kalign/) | [2.03](http://msa.sbc.su.se/downloads/kalign/current.tar.gz) | sequence alignment | download & compile | [16343337](https://www.ncbi.nlm.nih.gov/pubmed/16343337) |
| [T-Coffee](http://www.tcoffee.org/) | [11.00.8cbe486](http://www.tcoffee.org/Packages/Stable/Latest/T-COFFEE_distribution_Version_11.00.8cbe486.tar.gz) | alignment combining | conda | [10964570](https://www.ncbi.nlm.nih.gov/pubmed/10964570)
| [trimAl](http://trimal.cgenomics.org/) | [1.4.1](https://github.com/scapella/trimal/archive/v1.4.1.tar.gz) | alignment trimming | conda | [19505945](https://www.ncbi.nlm.nih.gov/pubmed/19505945)
||
| [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html) | [8.2.9](https://github.com/stamatak/standard-RAxML/archive/v8.2.9.tar.gz) | tree building | conda | [24451623](https://www.ncbi.nlm.nih.gov/pubmed/24451623) |
| [PhyML](http://www.atgc-montpellier.fr/phyml/) | [3.2.20160530](https://github.com/stephaneguindon/phyml/archive/v3.2.20160530.tar.gz) | tree building | download & compile | [20525638](https://www.ncbi.nlm.nih.gov/pubmed/20525638) |
| [FastTree](http://www.microbesonline.org/fasttree/) | [2.1.9](http://www.microbesonline.org/fasttree/FastTree) | tree building | conda | [20224823](https://www.ncbi.nlm.nih.gov/pubmed/20224823) |
||
| [MCL](http://micans.org/mcl/) | [14-137](http://micans.org/mcl/src/mcl-14-137.tar.gz) | clustering | conda | [22144159](https://www.ncbi.nlm.nih.gov/pubmed/22144159) |
| [Prodigal](http://prodigal.ornl.gov/) | [2.6.3](https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz) | gene calling | conda | [20211023](https://www.ncbi.nlm.nih.gov/pubmed/20211023) |

### Appendix: Availability of applications

#### Availability from conda channels

Multiple applications can be directly installed from conda channels. This not only simplifies the installation process but also guarantees the modularity of the entire pipeline. That is, the installed applications are only available from within the conda environment.

Available in [bioconda](https://anaconda.org/bioconda):
* mcl=14.137
* t_coffee=11.0.8
* muscle=3.8.1551
* raxml=8.2.9
* trimal=1.4.1
* diamond=0.8.24 (Note: DIAMOND in bioconda is newer than that in biocore)
* blast-legacy=2.2.22
* blast=2.2.31
* mafft=7.305
* orthofinder=1.1.2 (Note: OrthoFinder is only compatible with Python 2)
* prodigal=2.6.2 (Note: the default version is 2.60, which is older)
* fasttree=2.1.9 (Note: this install has double precision support, which is non-default but critical ([details](http://darlinglab.org/blog/2015/03/23/not-so-fast-fasttree.html)))

Available in [biocore](https://anaconda.org/biocore):
* prodigal=2.6.2
* blast-legacy=2.2.22
* blast-plus=2.2.31
* diamond=0.7.10

#### Availability from Ubuntu repositories

Ubuntu users may take advantage of the repositories. For example, **PhyML** and **Kalign**, which are not available from conda, can be installed by:
```
sudo apt-get install phyml kalign
```
Please note that the installation takes effect system-wide.

Available from the default Ubuntu 16.04 LTS repository (universe):
* fasttree (2.1.8)
* kalign (2.03+20110620)
* mafft (7.271)
* mcl (14-137)
* muscle (3.8.31)
* mysql-server (5.7.17)
* ncbi-blast+ (2.2.31)
* phyml (3.2.0)
* prodigal (2.6.2)
* raxml (8.2.4)
* t-coffee (11.00.8cbe486)

#### Cross-dependencies

* OrthoFinder requires MCL, FastTree, Blast+ and MAFFT
* Phylomizer requires Blast+, KAlign, MAFFT, MUSCLE, T-Coffee, and PhyML
* DarkHorse requires MySQL

#### Additional notes

The [ETE team](http://etetoolkit.org/) has released a conda package: [ete3_external_apps](https://anaconda.org/etetoolkit/ete3_external_apps), which wraps up multiple popular phylogenetics applications (including PhyML and Kalign). This is not required for Horizomer, but if you wish to install, do:
```
conda install -c etetoolkit ete3_external_apps
```

(WIP: A script will be created to automate the installation of most of these applications. However, a few of them are definitely not automatable due to license issues.)
