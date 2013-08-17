#!/usr/bin/perl

#fasta2popgen.pl
#script to calculate various pop gen metrics per position from aligned fasta entries

#Megan Supple
#last modified 14 June 2013

#usage fasta2popgen.pl [options] <in.fasta> <in.pops>
	#options are to specifiy the analysis 
		#fst (based on Weir 1996, Genetic Data Analysis II) 
		#3 level hierarchical fst (by population, then pheontype) (based on Weir 1996, Genetic Data Analysis II)
		#genotype by phenotype association (two tailed Fisher's exact test)
		#measures of selection (heterozygosity, pi, segregating sites, Tajima's D, D*, singletons, SNP type (private vs shared)

#requires BioPerl, FormatGenos.pm, PopGen.pm, PopStatsHierarchy.pm, PopStatsModified.pm

#produces per position values for each analysis (see PopGen modules for details)


use lib $ENV{PERL5LIB};
use strict;
use warnings;

use FormatGenos;
use PopGen;
use Getopt::Long;
use List::MoreUtils qw(firstidx);
use Data::Dumper;

my $usage = "Usage: fasta2popgen.pl [options] <in.fasta> <in.pops>
options:
	<in.fasta> 	an aligned fasta file
	<in.pops>	a space delimited text file with information on each sample \n\t\t\t (column 1: population, column 2: phenotype, column 3: name in input fasta) 
	-fst		calculate Fst (samples grouped by pheotype, ignoring population information)
	-fst3level	calculate 3 level fst (level 1: population, level 2: phenotype), returns fst based on pheontype
	-assoc		calculate genotype by phenotype association (samples grouped by pheotype, ignoring population information)
	-sel		calculate various measures used to test for selection (samples grouped by phenotype, ignoring population information)
";


my $calc_fst;
my $calc_fst3;
my $calc_gxp;
my $calc_sel;

GetOptions (   "fst" => \$calc_fst,
               "fst3level" => \$calc_fst3,
               "assoc" => \$calc_gxp,
               "sel" => \$calc_sel,
	    );


#read in command line argument
my ($infasta,$inpops)=@ARGV;
die "$usage" unless (@ARGV == 2);


########################################################################
###parse input files####################################################
########################################################################

#initialize hash for samples and population information
my %samples = ();

#open and read in input file with sample and population information
open(INPOPS, $inpops)||die "can't open input population file. $!\n";
my @input=<INPOPS>;
close INPOPS;
my $n_samples=@input;

#process each sample's population information
for (my $i=0; $i<@input; $i++)
	{
	  #split information
	  my @info=split(" ",$input[$i]);
	  #enter information into the sample hash
	  $samples{$info[2]}={'pop' => $info[0], 'pheno' => $info[1]};
	  #print "population= $samples{$info[2]}{'pop'}\n";
	  #print "phenotype= $samples{$info[2]}{'pheno'}\n"; 
	}

#parse input fasta
my ($seqs_p,$contig,$size)=FormatGenos::fasta2hash($infasta);


#########################################################################
###add sequence info to data structure###################################
#########################################################################

#initialize pointers to population and phenotype groupings
my @pheno_names=();             #array with phenotype names
my @phenos;
my @pheno_header=();
my @pop_names=();
my @pops;
my @pop_header=();
my $n_total=0;

#process each fasta entry if the sample exists in the list
while (my ($seqid, $seq) = each(%$seqs_p))
	{
	  #if the sample exist in the sample list
	  if (exists $samples{$seqid})
		{
		  #get population and phenotype information
		  my $pop=$samples{$seqid}{'pop'};
		  my $pheno=$samples{$seqid}{'pheno'};
		  if (!defined $pop){die "sample $seqid does not exist in population input file $inpops\n"}
	
		  #process the entry
		  print "processing $seqid\n";
		  my $indiv_p=FormatGenos::seq2geno($seqid,$seq);


		  #if grouping "population" by pheontype (for fst, association, and selection)
		  if ($calc_fst || $calc_gxp || $calc_sel)
	        	{	  
			  #create pheotype group (if it is new) and push into @phenos
			  if (!grep(/^$pheno$/,@pheno_names))
				{
				  push(@pheno_names, $pheno);
				  push(@phenos, Bio::PopGen::Population->new(-name => $pheno));
				}

			  #add individual to the phenotypic group by finding the index of the phenotype designation
			  $phenos[firstidx {$_ eq $pheno} @pheno_names]->add_Individual($$indiv_p);
			}
		  #if grouping first by popualation and then by phenotype (fst3level)
		  if ($calc_fst3)
		        {
			  #create population (if it is new) and push into @pops
			  if (!grep(/^$pop\_$pheno$/,@pop_names))
				{
				  push(@pop_names, "$pop\_$pheno");
				  push(@pops, Bio::PopGen::Population->new(-name => "$pop\_$pheno", -description => $pheno, -source => $pop));
				}
			  #add individual to the population group by finding the index of the phenotype designation
	                  $pops[firstidx {$_ eq "$pop\_$pheno"} @pop_names]->add_Individual($$indiv_p);
		        }
		}
	  #if sample is not in list, ignore it
	  print "sample $seqid is not in the sample list\n"; 
	}

#make headers
#if grouping "population" by pheontype (for fst, association, and selection)
if ($calc_fst || $calc_gxp || $calc_sel)
  	{
	  @pheno_header=("####################\n");
	  #iterate over phenotypes
	  for (my $i=0; $i<@pheno_names; $i++)
		{
		  #determine sample size for the phenotype and print to header
		  my $num_inds=$phenos[$i]->get_number_individuals;
		  push(@pheno_header,"#$pheno_names[$i]:n=$num_inds\n#");
		  #get a list of samples associated with the phenotype
		  my @inds=$phenos[$i]->get_Individuals();
		  foreach(@inds){push(@pheno_header, $_->unique_id);push(@pheno_header,","); }
		  push(@pheno_header, "\n");
		  $n_total+=$num_inds;
		}
	  push(@pheno_header,"####################\n");
	}
#if grouping by population, then by phenotype (for fst3level)
if ($calc_fst3)
	{
	  @pop_header=("####################\n");
          #iterate over phenotypes
          for (my $i=0; $i<@pop_names; $i++)
		{
		  #determine sample size for the population and print to header
                  my $num_inds=$pops[$i]->get_number_individuals;
                  push(@pop_header,"#$pop_names[$i]:n=$num_inds\n#");
                  #get a list of samples associated with the population
                  my @inds=$pops[$i]->get_Individuals();
                  foreach(@inds){push(@pop_header, $_->unique_id);push(@pop_header,","); }
                  push(@pop_header, "\n");
                  $n_total+=$num_inds;
                }
           push(@pop_header,"####################\n");
	}
	


#check if all individuals in input file have been entered
if ($n_total != $n_samples){die "The number of samples in $inpops does not equal the number of samples in $infasta";}

###############################################################################
###calculate popgen estimates##################################################
###############################################################################

#calculate Fst
if ($calc_fst)
	{
	  print "calculating Fst\n";
	  PopGen::Fst(\@pheno_names, \@phenos, \$contig, \$size, \@pheno_header);
	}

#calculate 3 level Fst
if ($calc_fst3)
        {
          print "calculating 3 level Fst\n";
          PopGen::Fst_3level(\@pop_names, \@pops, \$contig, \$size, \@pop_header);
        }

#calculate genotype by phenotype association
if ($calc_gxp)
	{
	  print "calculating association\n";
	  PopGen::association(\@pheno_names, \@phenos, \$contig, \$size, \@pheno_header);
	}

#calculate selection measures
if ($calc_sel)
	{
	  print "calculating measures of selection\n";
	  PopGen::selection(\@pheno_names, \@phenos, \$contig, \$size, \@pheno_header);
	}


print "done!\n";

