#!/bin/bash

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, The WGS-HGT Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# usage: automate the installation of DarkHorse2
# this script will install:
#   MySQL 5.7.18
#   NCBI nr, taxdump and prot.accession2taxid (up-to-date)
#   DarkHorse 2.0_rev06
# the whole system (including MySQL) will run in a local directory and won't
# affect anything outside

set -eu
appdir=$(readlink -m $1)

# download and compile MySQL source code
mkdir -p $appdir/mysql
wget https://dev.mysql.com/get/Downloads/MySQL-5.7/mysql-boost-5.7.18.tar.gz
tar xvf mysql-boost-5.7.18.tar.gz
cd mysql-5.7.18
cmake -D MYSQL_DATADIR=$appdir/mysql/data \
      -D SYSCONFDIR=$appdir/mysql/etc \
      -D CMAKE_INSTALL_PREFIX=$appdir/mysql \
      -D WITH_BOOST=./boost .
make
make install

# clean up installation files
cd ..
rm mysql-boost-5.7.18.tar.gz
rm -rf mysql-5.7.18

# find an unused port number
#   the default MySQL port is 3306
#   to avoid conflict with system-wide MySQL (if any), search from 3307
port=3307
while true
do
    netstat -tlnp | grep $port && ((port++)) || break
done

# create MySQL configuration file
mkdir -p $appdir/mysql/etc
echo "
[client]
port           = $port
socket         = $appdir/mysql/mysql.sock
[mysqld]
port           = $port
socket         = $appdir/mysql/mysql.sock
basedir        = $appdir/mysql
datadir        = $appdir/mysql/data
sql_mode       = NO_ENGINE_SUBSTITUTION
[mysqld_safe]
log-error      = $appdir/mysql/data/mysql.err
pid-file       = $appdir/mysql/mysql.pid
" > $appdir/mysql/etc/my.cnf

export PATH=$appdir/mysql/bin:$PATH

# initialize MySQL serivce
mysqld --initialize-insecure
sleep 60

# start MySQL serivce
mysqld &
sleep 60

# set up MySQL account
#   username: root, password: password
mysql -u root --skip-password -e "ALTER USER 'root'@'localhost' IDENTIFIED BY 'password';"

# download NCBI database
mkdir -p $appdir/ncbi
cd $appdir/ncbi
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
gunzip nr.gz
mkdir -p taxdump
tar xvf taxdump.tar.gz -C taxdump
rm taxdump.tar.gz
gunzip prot.accession2taxid.gz

# install required Perl modules
export PERL_MM_USE_DEFAULT=1
perl -MCPAN -e 'install App::cpanminus'
cpanm DBI DBD::mysql

# download DarkHorse 2-rev06
dhv=rev06
mkdir -p $appdir/darkhorse
cd $appdir/darkhorse
wget https://github.com/spodell/Darkhorse2/archive/$dhv.tar.gz
tar xvf $dhv.tar.gz
mv Darkhorse2-$dhv 2-$dhv
rm $dhv.tar.gz

# patch DarkHorse scripts
#   let DarkHorse connect MySQL via a custom socket
before='my $db_path = "DBI:$db_program:$db_name:$db_host:";'
after1='my $db_path = "DBI:$db_program:$db_name:$db_host;mysql_socket=".$config_params{db_socket};'
after2='my $db_path = "DBI:$db_program:$db_name:$db_host;mysql_socket=".$dh_config{db_socket};'
cd $appdir/darkhorse/2-$dhv
while read pl
do
    chmod +x $pl
    grep -q "%config_params" $pl \
         && sed -i "s/$before/$after1/g" $pl \
         || sed -i "s/$before/$after2/g" $pl
done <<< "$(find . -name '*.pl' -type f)"

# create DarkHorse configuration file
echo "
# DarkHorse2 configuration file for the WGS-HGT pipeline
# installation parameters
[program_directory]=$appdir/darkhorse/2-$dhv
[genbank_nr_fasta_path]=$appdir/ncbi/nr
[names_dmp_path]=$appdir/ncbi/taxdump/names.dmp
[nodes_dmp_path]=$appdir/ncbi/taxdump/nodes.dmp
[protid_taxid_dmp_path]=$appdir/ncbi/prot.accession2taxid
# database access parameters
[db_name]=dh2nr
[db_program]=mysql
[db_host]=localhost
[db_user]=root
[db_user_password]=password
[db_socket]=$appdir/mysql/mysql.sock
[max_lines_per_packet]=4000
[min_lineage_terms]=2
[min_align_coverage]=0.7
" > $appdir/darkhorse/config.txt

# create DarkHorse database
mysql -u root -ppassword -e "CREATE DATABASE dh2nr;"

# run DarkHorse installation script
perl install_darkhorse2.pl -c $appdir/darkhorse/config.txt

# stop MySQL service
mysqladmin -u root -ppassword shutdown

# finish up
echo "
Summary:
MySQL 5.7.18 has been installed at: $appdir/mysql.
  user: root
  password: password
  port: $port
  socket: $appdir/mysql/mysql.sock
DarkHorse 2-$dhv has been installed at: $appdir/darkhorse.
NCBI nr database has been downloaded at: $appdir/ncbi.
MySQL database dh2nr has been created.

To use DarkHorse, start MySQL first:
  export PATH=$appdir/mysql/bin:\$PATH
  mysqld &

When completed, shut down MySQL:
  mysqladmin -u root -ppassword shutdown

"
echo "Task completed."
