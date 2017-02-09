### Installation (Linux & Mac OS)

#### Recommended installation order

Have [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://anaconda.org/) installed in the system. The Python 3 version is recommended.

Create conda environment and install required libraries:
```
conda create -n wgshgt python=3.5
conda install click biopython
conda install -c biocore scikit-bio
```

Install third-party applications using conda:
```
conda install -c bioconda blast blast-legacy diamond fasttree mafft mcl muscle orthofinder prodigal=2.6.2 raxml t_coffee trimal
```

Install more third-party applications:

If you are a Ubuntu user, do:
```
sudo apt install phyml kalign
```
Otherwise, download and install these programs manually. See the application list below.

Tip: You can used this command to install multiple programs automatically, including FastTree, RAxML, Phyml, Trimal, Pmodeltest, T-Coffee, M-Coffee, Kalign, Prank, Probcons, Muscle, ClustalOmega, Dialign-tx, Mafft, Consel, SLR and Codeml:
```
conda install -c etetoolkit ete3_external_apps
```


#### List of third-party applications

| Name | Tested Version | Purpose | Required | Auto-install | PMID |
| --- | --- | --- | --- | --- | --- |
| T-REX | 3.4 | HGT detection | yes | yes | |
| RANGER-DTL-U | [1.0](http://compbio.mit.edu/ranger-dtl/ranger-dtl-linux.tar.gz) | HGT detection | yes | yes | |
| PhyloNet | [3.5.6](http://bioinfo.cs.rice.edu/sites/bioinfo.cs.rice.edu/files/kcfinder/files/PhyloNet_3.5.6.jar) | HGT detection | yes | yes | |
| Jane | [4](https://www.cs.hmc.edu/~hadas/jane/index.html) | HGT detection | yes | yes | |
| TREE-PUZZLE | [5.3.rc16](http://www.tree-puzzle.de/tree-puzzle-5.3.rc16-linux.tar.gz) | HGT detection | yes | yes | |
| CONSEL | [1.20](http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/pub/cnsls020.tgz) | HGT detection | yes | yes | |
| DarkHorse | [1.5 rev170](http://darkhorse.ucsd.edu/DarkHorse-1.5_rev170.tar.gz) | HGT detection | yes | yes | [17274820](https://www.ncbi.nlm.nih.gov/pubmed/17274820) |
| HGTector | [0.2.1](https://github.com/DittmarLab/HGTector/archive/wgshgt.zip) | HGT detection | yes | yes | [25159222](https://www.ncbi.nlm.nih.gov/pubmed/25159222) |
| EGID | [1.0](http://www5.esu.edu/cpsc/bioinfo/software/EGID/EGID_1.0.tar.gz) | HGT detection | yes | yes | [22355228](https://www.ncbi.nlm.nih.gov/pubmed/22355228) |
| GeneMarkS | [4.30](http://exon.gatech.edu/GeneMark/license_download.cgi) | HGT detection | yes | no | [9461475](https://www.ncbi.nlm.nih.gov/pubmed/9461475) |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) | [1.1.4](https://github.com/davidemms/OrthoFinder/releases/download/1.1.4/OrthoFinder-1.1.4.tar.gz) | orthology identification | yes | yes | [26243257](https://www.ncbi.nlm.nih.gov/pubmed/26243257) |
| [Phylomizer](https://github.com/Gabaldonlab/phylomizer) | [9/12/2016](https://github.com/Gabaldonlab/phylomizer/commit/e427a04b3d62bbac4d760fef975f6bdf5aeed44a) | gene family tree building | yes | yes | NA |
| [PhyloPhlAn](https://github.com/davidemms/OrthoFinder) | [1/25/2017](https://bitbucket.org/nsegata/phylophlan/commits/2c0e61ad820b8ff732837e98f7843afbb7ec1cda) | genome tree building | yes | yes | [23942190](https://www.ncbi.nlm.nih.gov/pubmed/23942190) |
||
| [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) | [2.5.0](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/ncbi-blast-2.5.0+-src.tar.gz) | sequence similarity search | yes | yes | [2231712](https://www.ncbi.nlm.nih.gov/pubmed/2231712) |
| [DIAMOND](https://ab.inf.uni-tuebingen.de/software/diamond) | [0.8.28](https://github.com/bbuchfink/diamond/archive/v0.8.28.tar.gz) | sequence similarity search | yes | yes | [25402007](https://www.ncbi.nlm.nih.gov/pubmed/25402007) |
||
| [MAFFT](http://mafft.cbrc.jp/alignment/software/) | [7.305](http://mafft.cbrc.jp/alignment/software/mafft-7.305-with-extensions-src.tgz)  | sequence alignment | yes | yes | [12136088](https://www.ncbi.nlm.nih.gov/pubmed/12136088) |
| [MUSCLE](http://drive5.com/muscle/) | [3.8.31](http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz) | sequence alignment | yes | yes | [15034147](https://www.ncbi.nlm.nih.gov/pubmed/15034147) |
| [KAlign](http://www.ebi.ac.uk/Tools/msa/kalign/) | [2.03](http://msa.sbc.su.se/downloads/kalign/current.tar.gz) | sequence alignment | no | yes | [16343337](https://www.ncbi.nlm.nih.gov/pubmed/16343337) |
| [T-Coffee](http://www.tcoffee.org/) | [11.00.8cbe486](http://www.tcoffee.org/Packages/Stable/Latest/T-COFFEE_distribution_Version_11.00.8cbe486.tar.gz) | alignment combining | yes | yes | [10964570](https://www.ncbi.nlm.nih.gov/pubmed/10964570)
||
| [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html) | [8.2.9](https://github.com/stamatak/standard-RAxML/archive/v8.2.9.tar.gz) | tree building | yes | yes | [24451623](https://www.ncbi.nlm.nih.gov/pubmed/24451623) |
| [PhyML](http://www.atgc-montpellier.fr/phyml/) | [3.2.20160530](https://github.com/stephaneguindon/phyml/archive/v3.2.20160530.tar.gz) | tree building | yes | yes | [20525638](https://www.ncbi.nlm.nih.gov/pubmed/20525638) |
| [FastTree](http://www.microbesonline.org/fasttree/) | [2.1.9](http://www.microbesonline.org/fasttree/FastTree) | tree building | yes | yes | [20224823](https://www.ncbi.nlm.nih.gov/pubmed/20224823) |
||
| [MCL](http://micans.org/mcl/) | [14-137](http://micans.org/mcl/src/mcl-14-137.tar.gz) | clustering | no | yes | [22144159](https://www.ncbi.nlm.nih.gov/pubmed/22144159) |
| [Prodigal](http://prodigal.ornl.gov/) | [2.6.3](https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz) | gene calling | yes | yes | [20211023](https://www.ncbi.nlm.nih.gov/pubmed/20211023) |

#### Conda install
The following programs can be installed directly from conda:

Available in [bioconda](https://anaconda.org/bioconda):
* mcl=14.137
* t_coffee=11.0.8
* muscle=3.8.1551
* raxml=8.2.9
* trimal=1.4.1
* diamond=0.8.24 (Note: diamond in bioconda is newer than that in biocore)
* blast-legacy=2.2.22
* blast=2.2.31
* mafft=7.305
* orthofinder=1.1.2
* prodigal=2.6.2 (Note: the default version is 2.60, which is older)
* fasttree=2.1.9

Available in [biocore](https://anaconda.org/biocore):
* prodigal=2.6.2
* blast-legacy=2.2.22
* blast-plus=2.2.31
* diamond=0.7.10

#### Cross-dependencies
* OrthoFinder requires MCL, FastTree, Blast+ and MAFFT
* Phylomizer requires Blast Legacy, KAlign, MAFFT, MUSCLE, T-Coffee, PhyML
