#!/bin/bash

# this script will install:
#   MySQL 5.7.17
#   NCBI nr, taxdump and prot.accession2taxid
#   DarkHorse 2.0_rev05

appdir=$(readlink -m $1)

# download and compile MySQL source code
mkdir -p $appdir/mysql
wget https://dev.mysql.com/get/Downloads/MySQL-5.7/mysql-boost-5.7.17.tar.gz
tar xvf mysql-boost-5.7.17.tar.gz
cd mysql-5.7.17
cmake -D MYSQL_DATADIR=$appdir/mysql/data \
      -D SYSCONFDIR=$appdir/mysql/etc \
      -D CMAKE_INSTALL_PREFIX=$appdir/mysql \
      -D WITH_BOOST=./boost .
make
make install

# clean up installation files
cd ..
rm mysql-boost-5.7.17.tar.gz
rm -rf mysql-5.7.17

export PATH=$appdir/mysql/bin:$PATH

# initialize MySQL serivce
#   username: root, password: (no password)
mysqld --initialize-insecure \
       --datadir=$appdir/mysql/data

# start MySQL serivce
#   3306 is the default MySQL port.
#   to avoid conflict with the system-wide MySQL, 3307 is used here.
mysqld --basedir=$appdir/mysql \
       --datadir=$appdir/mysql/data \
       --log-error=$appdir/mysql/data/mysql.err \
       --pid-file=$appdir/mysql/mysql.pid \
       --socket=$appdir/mysql/mysql.sock \
       --port=3307 &

# One may connect to the MySQL server with:
#  mysql -S $appdir/mysql/mysql.sock -u root --skip-password

# create DarkHorse database
mysql -u root --skip-password -e "CREATE DATABASE Darkhorse2_01;"

# stop MySQL service
mysqladmin -S $appdir/mysql/mysql.sock -u root --skip-password

# download NCBI database
mkdir -p $appdir/ncbi
cd $appdir/ncbi
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
gunzip prot.accession2taxid.gz

# install required Perl modules
export PERL_MM_USE_DEFAULT=1
perl -MCPAN -e 'install App::cpanminus'
cpanm DBI DBD::mysql

mkdir -p $appdir/darkhorse
cd $appdir/darkhorse

# create DarkHorse configuration file
echo "
# DarkHorse2 configuration file for the WGS-HGT pipeline

# installation parameters
[program_directory]=$appdir/darkhorse/2-rev_05
[genbank_nr_fasta_path]=$appdir/ncbi/nr
[names_dmp_path]=$appdir/ncbi/taxdump/names.dmp
[nodes_dmp_path]=$appdir/ncbi/taxdump/nodes.dmp
[protid_taxid_dmp_path]=$appdir/ncbi/prot.accession2taxid

# database access parameters
[db_name]=Darkhorse2_01
[db_program]=mysql
[db_host]=localhost
[db_user]=root
[db_user_password]=
[max_lines_per_packet]=32000
[min_lineage_terms]=2
[min_align_coverage]=0.7
" > $appdir/darkhorse/config.txt

# download and install DarkHorse2
wget https://github.com/spodell/Darkhorse2/archive/rev_05.tar.gz
tar xvf rev_05.tar.gz
mv Darkhorse2-rev_05 2-rev_05
cd 2-rev_05
perl install_darkhorse2.pl -c $appdir/darkhorse/config.txt
cd ..
rm rev_05.tar.gz

echo done
