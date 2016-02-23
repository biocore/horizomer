#!/usr/bin/perl
# generate_dh_keywords.pl
# Sheila Podell
# February 15, 2016

# Based on recursive SQL queries to NCBI taxdump names and nodes tables,
# creates DarkHorse self-exclusion lists 
# reports at genus, species, strain and primary tax_id taxonomic granularity

# Takes as input a tab-delimited list of genome sequences with the following fields:
	# genome name
	# ncbi_tax_id_number 
# lines starting with "#" are ignored
# Returns warning and skips entry if tax_id not found in NCBI taxdump database
# writes tab-delimited ouptut with genome name, taxid, # hits for primary, strain, species, genus tax ids

# results options:	
	# -primary    create primary exclusion list file for each genome 
	# -strain     create strain exclusion list file for each genome 
	# -species    create species exclusion list file for each genome
	# -genus      create genus exclusion list file for each genome
	
use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::Spec;
use DBI;

# get command line arguments
	my ($db, $datafile, $sth, $sql, $config_filename, $printline, $debug, $primary, 
		$strain, $species, $genus, $show_help, $result_dir);
	my $USAGE;
	my $message = qq(
  Usage: $0 -i source_datafile -c config filename 
  optional parameters:
  	-h           print help message and exit
  	-primary     output primary tax-id exclusion list file
  	-strain      output strain tax-id exclusion list file
  	-species     output species tax-id exclusion list file (includes strains)
  	-genus       output genus tax-id exclusion list file
  );
	GetOptions( 
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
			    "d=s" => \$debug,
			    'primary' => \$primary, 
			    'strain' => \$strain,
			    'species' => \$species,
			    'genus' => \$genus,
			    'help' => \$show_help,
			    'h' => \$show_help,
				);	     
				
# check arguments
	if($USAGE || !$datafile || $show_help)
	{
		print STDERR "$message\n";
		exit(0);
	} 	
	unless (defined $datafile && -s $datafile)
	{
		print STDERR "\n  Couldn't find source datafile.";
		print STDERR "$message\n";
		exit(0);
	}	
	unless (defined $config_filename && -s $config_filename)
	{
		print STDERR "\n  Couldn't find MySQL config file.";
		print STDERR "$message\n";
		exit(0);
	}

#global variables
	my $warning_count = 0;
	my $error_count = 0;
	my @exclude_list_template = ("cloning",
"vector",
"plasmid",
"cosmid",
"expression",
"environmental",
"synthetic",
"construct",
"contaminant",
"unidentified",
"unknown",
"untyped",
"unspecified",
"clone",); 

# get SQL database parameters from config file
	my %config_params = &read_config_file($config_filename);
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";
	my $min_lineage_terms = $config_params{min_lineage_terms} || "not found";
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";
	my @row = ();
	my $element;
	print STDERR "  Connecting to $db_name database . . .\n";

#set up a working directory, move there
	$result_dir = "$$"."_results";
	mkdir "$result_dir", 0755 or warn "Cannot make results directory $result_dir\n: $!";
	my $cwd = getcwd();
		
# open a logfile, put STDERR comments there too.
	my $logfilename = "$$"."_logfile";
	my $date = `date`;
	chomp $date;
	open (LOGFILE, ">$result_dir/$logfilename") or die "can't open log file $logfilename\n $!\n";
	print LOGFILE qq(begin: $date
$0
input file = $datafile
config file = $config_filename 
);

# open a file to hold summary of results
	my $tabname = "$$"."_results_smry";
	open (TABFILE, ">$result_dir/$tabname") or die "can't open tab summary file $tabname\n $!\n";
	my @col_headers = ("genome_name", "primary_tax_id", "unix_name", "num_genus_ids", "num_species_ids", "num_strain_ids",);
	my $headers = join "\t", @col_headers;	
	print TABFILE "$headers\n";
	
# get list of input names, create genome objects
	my @genome_list = ();		# original user entry order
	my %genome_objects = ();	# key = user-entered genome name, value = object ref
	my %uniq_primary_ids =();
	my $current = "";
    open (INPUT, "$datafile")  or die "can't input file $datafile\n$!\n";
    while (<INPUT>)
    {
   		chomp;
   		next if $_ =~ /^\s*$/;	# ignore blank lines
   		next if $_ =~ /^#/;	    # ignore comments
   		next if $_ =~ /^genome/;	# ignore header line
   		my @tmp = split "\t", $_; 		
   		unless (scalar @tmp == 2)
   		{
   			print STDERR "Input data error, line $.:$_\n";
   			print LOGFILE "Input data error, line $.:\n$_\n";
   			print STDERR "Expected format: species name<tab>ncbi_tax_id \n";
   			$error_count++;
   			next;
   		}  		
   	# warn about duplicate names
   		my $genome_name = $tmp[0];
   		my $primary_tax_id =  $tmp[1];
   		if (exists $genome_objects{$genome_name})
		{
			print STDERR "  WARNING: input contains duplicate genome name: $genome_name\n";
			print LOGFILE "  WARNING: input contains duplicate genome name: $genome_name\n";
			$warning_count++;
			next;
		}
		
	# don't try to process genomes without ncbi tax id
		unless (defined $primary_tax_id && ($primary_tax_id =~/\d+/) && $primary_tax_id >0) 
		{
			print STDERR "  ERROR: skipping $genome_name: no primary tax_id found\n";
			print LOGFILE "  ERROR: skipping $genome_name: no primary tax_id found\n";
			$primary_tax_id = 0;
			$error_count++;
			next;
		}
		
	# check to be sure the primary tax_id is in the database
		my $tax_id_check = &find_tax_id_in_db_version($primary_tax_id);
		unless (defined $tax_id_check && $tax_id_check > 0)
		{
			print STDERR "  ERROR: skipping $genome_name: primary tax_id not in database $db_name\n";
			print LOGFILE "  ERROR: skipping $genome_name: primary tax_id not in database $db_name\n";
			$primary_tax_id = 0;
			$error_count++;
			next;
		}
		
	# check for duplication
		if (defined $primary_tax_id && $primary_tax_id > 0)
		{
			$uniq_primary_ids{$primary_tax_id}++;
			if ($uniq_primary_ids{$primary_tax_id} > 1)
			{
				#print STDERR "  WARNING: non-unique primary tax id $primary_tax_id, $genome_name\n";
				print LOGFILE "  WARNING: non-unique primary tax id $primary_tax_id, $genome_name\n";
				$warning_count++;
			}
		}
	
	# create object to hold genome attributes
		my @args = ($genome_name, $primary_tax_id);
		my $current = &new_genome("genome", \@args);
		if (defined $current) 
		{
			$genome_objects{$genome_name} = $current;
			push @genome_list, $genome_name;
		}	
	}
  close INPUT;
  
# Count number of genomes, give feedback
	my $num_genomes = scalar keys %genome_objects; 
	unless ($num_genomes >0)
	{
		print STDERR "ERROR: no valid data to process\n";
		print LOGFILE "ERROR: no valid data to process\n";
		exit(0);
	}	  
 
# get taxonomic siblings from SQL databases 
	my $num_analyzed = 0;
	foreach my $genome_name (@genome_list)
	{
		my ($num_species, $num_strains);
		my $genome_obj = $genome_objects{$genome_name};
		my $primary_tax_id = $genome_obj->{primary_tax_id};
		my $current_id;
		my $current_name = $genome_name;
		my $rank = "";
		my $genus_name = "";
		my $counter = 0;  # variable to prevent run-away recursion
		
		$genome_obj->{primary_tax_id} = $primary_tax_id;		
	
	# get query-specific data by traversing linked lists in ncbi nodes table 
		$counter = 0;
		$current_id = $primary_tax_id;
		
		while ($counter < 4)
		{				
			$counter++;	
			my $result_ref = &get_tax_id_rank($current_id);				
			my @result = @$result_ref;
			
			my $new_rank =  $result[0];
			my $new_name =  $result[1];
			my $new_name_class =  $result[2];										
			if (defined $new_rank && length $new_name > 3
				&& defined $new_rank && length $new_rank >3)
			{
				$rank = $new_rank;
				$current_name = $new_name;
			}
			else 			
			{
				print LOGFILE "  ERROR: undefined name/rank, input #$num_analyzed, $genome_name\n";
				last;
			}
			
			if ($rank eq "no rank" || $rank eq "subspecies") 
			{
				unshift @{$genome_obj->{subspecies_list}}, $current_id;
				$genome_obj->{strain_name} = $current_name;
				$genome_obj->{strain_id} = $current_id;	
			}
			if ($rank eq "species")
			{
				unshift @{$genome_obj->{species_list}}, $current_id; 
				$genome_obj->{species_name} = $current_name;
				 $genome_obj->{species_id} = $current_id;	
			}					
			if ($rank eq "species group")
			{
				my @tmp = split " ", $current_name;
				if ($current_name =~ /^Candidatus\s(.*)/)
				{
					shift @tmp;
				}
				$genus_name = $tmp[0];
			}
			if ($rank eq "genus")
			{		
				unshift @{$genome_obj->{genus_list}}, $current_id;
				if ($current_name =~ /^Candidatus\s(.*)/)
				{
					$current_name = $1;
				}
				my @tmp = split " ", $current_name;
				$genus_name = $tmp[0];
		
				# get the taxonomy id and the name
					$genome_obj->{genus_name} = $genus_name;
					unshift @{$genome_obj->{genus_list}}, $genus_name;										
			}
										
			last if ($rank eq "genus");
			my $next_id = &get_parent_tax_id($current_id);
			$current_id = $next_id;				
		}
		
		if ($num_analyzed == 0)
		{
			print STDERR "  Searching taxonomy id's for $num_genomes genomes . . .\n";
		}
		$num_analyzed++;
		if ($num_analyzed % 50 == 0)
		{
				print STDERR "  $num_analyzed genomes processed . . .\n";
		}
		&write_output($genome_name);
		print LOGFILE "   completed $genome_name\n";						
	}

$date = `date`;
print STDERR "  Sucessfully processed $num_analyzed genomes\n";
print STDERR "    $warning_count WARNINGS\n";
print STDERR "    $error_count ERRORS\n";
print STDERR "  See logfile $logfilename for details.\n";

print LOGFILE "  Sucessfully processed $num_analyzed genomes\nend: $date\n";
print LOGFILE "    $warning_count WARNINGS\n";
print LOGFILE "    $error_count ERRORS\n";

close LOGFILE;
close TABFILE;

################################################
# SUBROUTINES
################################################
sub read_config_file
{
    my ($filename) = @_;
    open (INPUT, $filename)  or die "can't open DarkHorse config file $filename\n$!\n";
	my %dh_config = ();
	while (<INPUT>)
	{
		next if ($_ =~ /^\s+$/);
		next if ($_ =~ /^#/);
		if ($_ =~ /\[(.+)\]=\s*(.+)/)
		{
			$dh_config{$1}=$2;		
		}			
	}
	close INPUT;			
    return %dh_config;
}

sub new_genome
{
 my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;  

# Constructor attributes 
	$self->{genome_name} = "$properties[0]";
	$self->{primary_tax_id} = "$properties[1]";
	
	my $unix_name = $self->{genome_name};
# get rid of spaces and inconvenient characters
	if ($unix_name =~ /(.+)\s+$/)
		{
			$unix_name = $1;
		}
	if ($unix_name =~ /^\s*(.+)$/)
		{
			$unix_name = $1;
		}
	$unix_name =~ tr/ /_/s;
	$unix_name =~ tr/-/_/;
	$unix_name =~ s/\[//g;
	$unix_name =~ s/\]//g;
	$unix_name =~ s/'//g;
	$unix_name =~ s/"//g;
	$unix_name =~ s/\(//g;
	$unix_name =~ s/\)//g;
	$unix_name =~ s/\.//g;
	$unix_name =~ s/,//g;
	$unix_name =~ s/://g;
	$unix_name =~ s/;//g;
	$unix_name =~ s/\///g;
	$unix_name =~ s/\\//g;
	$unix_name =~ s/#//g;
	$unix_name =~ s/&//g;
	$unix_name =~ s/\?//g;
		
	$self->{unix_name} = $unix_name;
	$self->{genus_name} = "";
	
	my @primary_list =  ("$properties[1]", @exclude_list_template);
	my @strain_list = ("$properties[1]", @exclude_list_template);
	my @species_list = ("$properties[1]",@exclude_list_template);
	my @genus_list = ("$properties[1]",@exclude_list_template);

	$self->{primary_list} = \@primary_list;
	$self->{strain_list} = \@strain_list;
	$self->{species_list} = \@species_list;
	$self->{genus_list} = \@genus_list;

   return($self)
}

sub get_tax_id_rank
{			
	my ($tax_id) = @_;
	my $current_rank;
	my $current_name;
	$sql = qq(select rank, name_txt, name_class
from names, nodes
where names.tax_id = nodes.tax_id
and names.tax_id = $tax_id;);
	
	if ($debug) {print STDERR "\n$sql\n";}		
		$sth = $dbh->prepare($sql)
			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
		$sth->execute()
			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
		@row = ();
		while (@row = $sth->fetchrow_array)
		{
			next if $row[0] =~ /name_txt/;	#don't include SQL header
			next if $row[2] =~ /type\smaterial/;
			$current_rank = "$row[0]";
			$current_name = "$row[1]";
			last if ("$row[2]" eq "scientific name");
		}
		my @out_array = ($current_rank,$current_name);
		my $out_ref = \@out_array;
		
		unless (defined $current_rank && length $current_rank > 3
				&& defined $current_name && length $current_name >3)
		{
				print STDERR "ERROR: undefined name/rank for $tax_id\n";
				print LOGFILE "ERROR: undefined name/rank for $tax_id\n";
				print STDERR "$sql\n";
				print LOGFILE "$sql\n";
				$error_count++;
				#exit(0);
		}
						
		return $out_ref;	 
	}
	
sub find_tax_id_in_db_version
{			
	my ($input_tax_id) = @_;
	my $found_id; 
	$sql = "select tax_id from names where tax_id = $input_tax_id;";
	if ($debug) {print STDERR "\n$sql\n";}		
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

	@row = ();
	my $id_count = 0;
	while (@row = $sth->fetchrow_array)
	{
		unless (defined $row[0])
		{
			return(0);
		}
		next if $row[0] =~ /tax_id/;	#don't include SQL header
		$found_id = $row[0];				
	}	
	
	return $found_id;	 
}

sub get_parent_tax_id
{			
	my ($input_tax_id) = @_;
	my $parent = "not_found"; 
	$sql = "select parent_tax_id from nodes where tax_id = $input_tax_id;";
	if ($debug) {print STDERR "\n$sql\n";}		
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

	@row = ();
	my $id_count = 0;
	while (@row = $sth->fetchrow_array)
	{
		next if $row[0] =~ /tax_id/;	#don't include SQL header
		$parent = "$row[0]";
	}	
	return $parent;	 
}

sub get_siblings
{		

		my ($parent_id) = @_;
		my @sibs_list = ();
		
		unless (defined $parent_id && $parent_id > 0)
		{
			return undef;
		}
		$sql = qq(select distinct nodes.tax_id 
from names, nodes
where names.tax_id = nodes.tax_id 
and parent_tax_id = $parent_id;);
	if ($debug) {print STDERR "\n$sql\n";}	
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

	@row = ();
	my $id_count = 0;
	while (@row = $sth->fetchrow_array)
	{ 
		next if $row[0] =~ /tax_id/;	#don't include SQL header
		push @sibs_list, $row[0];
	}			
		
	return \@sibs_list;
}

sub write_output
{
	my ($key) = @_;
	
	my ($genus_list, $species_list, $strain_list);
	my $genome_name = $genome_objects{$key}->{genome_name};
	my $primary_tax_id = $genome_objects{$key}->{primary_tax_id};
	my $unix_name = $genome_objects{$key}->{unix_name};		
	my @genus_array = @{$genome_objects{$key}->{genus_list}};
	my @species_array = @{$genome_objects{$key}->{species_list}};
	my @strain_array = @{$genome_objects{$key}->{strain_list}};
	my @primary_array = @{$genome_objects{$key}->{primary_list}};
	
  # get ids for taxonomically related strains (subspecies or no rank with same parent as query species) 
		
		my $strain_sibs_ref = &get_siblings($genome_objects{$key}->{species_id});
		my @strain_sibs;
		if (defined $strain_sibs_ref)
		{ 
			@strain_sibs= @$strain_sibs_ref;
			push @strain_array, @strain_sibs;
			my $num_sibs = scalar @strain_array;
			if ($num_sibs > 1000)
				{
					#print STDERR "  WARNING: found $num_sibs strain-sibling hits for $genome_name\n";
					print LOGFILE "  WARNING: found $num_sibs strain-sibling hits for $genome_name\n";
					$warning_count++;
				}
			if ($debug)
			{print STDERR "found strain sibs for $key\n";}
		}	

	# get ids for taxonomically related species (same parent genous as query genome)
		my $species_sibs_ref = &get_siblings($genome_objects{$key}->{genome_id});
		my @species_sibs;
		if (defined $species_sibs_ref)
		{ 
			@species_sibs = @$species_sibs_ref;
			push @species_array, @species_sibs;
			my $num_sibs = scalar @species_array;
			if ($num_sibs > 1000)
			{
				#print STDERR "  WARNING: found $num_sibs species-sibling hits for $genome_name\n";
				print LOGFILE "  WARNING: found $num_sibs species-sibling hits for $genome_name\n";
				$warning_count++;
			}
		}		
		my $num_default = scalar @exclude_list_template;
		my $num_species = scalar @species_array - $num_default;
		my $num_strain = scalar @strain_array - $num_default;
		my $num_genus = scalar @genus_array - $num_default;									
		my @attributes = ($genome_name, $primary_tax_id, $unix_name, 
			$num_genus, $num_species, $num_strain, ); 		
		
	# print tab sumry output
		my $printstring = join "\t", @attributes;				
		print TABFILE "$printstring\n";
		
	# print individual list files 
		if ($primary)
		{
			unless (-e "$unix_name")
			{
				mkdir "$result_dir/$unix_name";
			}
			chdir "$unix_name";
			open (LISTFILE, ">$result_dir/$unix_name/primary_exclude_list") or die "can't open primary_exclude_list for $unix_name \n $!\n";			
			foreach my $next (@primary_array)
			{
				print LISTFILE "$next\n";
			}			
			close LISTFILE;		
		}
		 
		if ($strain)
		{
			unless (-e "$unix_name")
			{
				mkdir "$result_dir/$unix_name";
			}
			chdir "$unix_name";
			open (LISTFILE, ">$result_dir/$unix_name/strain_exclude_list") or die "can't open strain_exclude_list for $unix_name \n $!\n";			
			my %uniq_ids = ();
			foreach my $next (@strain_array)
			{
				$uniq_ids{$next}++;
				if ($uniq_ids{$next} <2)
				{
					print LISTFILE "$next\n";
				}
			}			
			close LISTFILE;
		
		}
		if ($species)
		{
			unless (-e "$unix_name")
			{
				mkdir "$result_dir/$unix_name";
			}
			chdir "$unix_name";
			open (LISTFILE, ">$result_dir/$unix_name/species_exclude_list") or die "can't open species_exclude_list for $unix_name \n $!\n";			
			my %uniq_ids = ();
			
		# append tax id numbers for strain sibs to species sibs 	
			@species_array = (@species_array,@strain_array);
			foreach my $next (@species_array)
			{
				$uniq_ids{$next}++;
				if ($uniq_ids{$next} <2)
				{
					print LISTFILE "$next\n";
				}
			}			
			close LISTFILE;		
		}
		if ($genus)
		{
			unless (-e "$unix_name")
			{
				mkdir "$result_dir/$unix_name";
			}
			chdir "$unix_name";
			open (LISTFILE, ">$result_dir/$unix_name/genus_exclude_list") or die "can't open genus_exclude_list for $unix_name \n $!\n";			
			my %uniq_ids = ();
		# put genus name in front of genus array
			my $genus_name =  $genome_objects{$key}->{genus_name};
			if (defined $genus_name && length $genus_name > 2)
			{
				unshift @genus_array, $genus_name;
			}
			
			foreach my $next (@genus_array)
			{
				$uniq_ids{$next}++;
				if ($uniq_ids{$next} <2)
				{
					print LISTFILE "$next\n";
				}
			}			
			close LISTFILE;		
		}
	}


__END__

	
############
# NOTES
############

mysql> select distinct name_class from names;
+---------------------+
| name_class          |
+---------------------+
| synonym             |
| scientific name     |
| in-part             |
| blast name          |
| genbank common name |
| authority           |
| misspelling         |
| type material       |
| equivalent name     |
| includes            |
| genbank synonym     |
| common name         |
| misnomer            |
| genbank acronym     |
| anamorph            |
| genbank anamorph    |
| teleomorph          |
| acronym             |
+---------------------+
18 rows in set (1.87 sec)

mysql> select distinct rank from nodes;
+------------------+
| rank             |
+------------------+
| no rank          |
| superkingdom     |
| genus            |
| species          |
| order            |
| family           |
| subspecies       |
| subfamily        |
| tribe            |
| phylum           |
| class            |
| forma            |
| suborder         |
| superclass       |
| subclass         |
| varietas         |
| kingdom          |
| subphylum        |
| superfamily      |
| infraorder       |
| infraclass       |
| superorder       |
| subgenus         |
| parvorder        |
| superphylum      |
| species group    |
| species subgroup |
| subtribe         |
| subkingdom       |
+------------------+
29 rows in set (0.68 sec)


select rank, name_txt, name_class
	from names, nodes
	where names.tax_id = nodes.tax_id
	and names.tax_id = 329726;
	
+---------+--------------------------------+-----------------+
| rank    | name_txt                       | name_class      |
+---------+--------------------------------+-----------------+
| no rank | Acaryochloris marina MBIC11017 | scientific name |
+---------+--------------------------------+-----------------+

my %bad_genus_names =(
"alpha" => 0,
"archaeon" => 0,
"bacterium " => 0,
"beta" => 0,
"butyrate" => 0,
"butyrate" => 0,
"candidate" => 0,
"Candidatus" => 0,
"delta" => 0,
"division" => 0,
"endosymbiont" => 0,
"epsilon" => 0,
"gamma" => 0,
"haloalkaliphilic" => 0,
"haloalkaliphilic" => 0,
"halophilic" => 0,
"marine" => 0,
"nanoarchaeote" => 0,
"primary" => 0,
"secondary" => 0,
"secondary" => 0,
"synthetic" => 0,
"unclassified" => 0,
"uncultured" => 0,
"unidentified" => 0,
"unidentified" => 0,
);
	