#!/usr/bin/perl

#fasta_filter.pl
#script to filter multi fasta for missing data (N)

#M Supple
#last modified 9 July 2015

#usage fasta_filter.pl <in.fasta> <threshold>
	#threshold is the % Ns that causes a site to be removed
		#ie if only want to remove sites that are all Ns, the threshold is 100

#requires BioPerl
	#http://www.bioperl.org/wiki/HOWTO:Multiple_Sequence_Alignment

#produces a fasta with sites with too much missing data removed


use lib $ENV{PERL5LIB};
use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

use Bio::AlignIO;
use Bio::SimpleAlign;
 

#read in command line arguments
my ($infasta,$thresh)=@ARGV;

# Use Bio::AlignIO to read in the alignment
my $in = Bio::AlignIO->new(-file => $infasta, -format => 'fasta');
my $aln = $in->next_aln();


#basic pre stats
my $numsamples=$aln->num_sequences;
print "number of samples is originally $numsamples\n";
my $seqlen=$aln->length;
print "sequence length is originally $seqlen\n";


#count Ns in each column and record column numbers with %Ns >= threshold
my $ncount;

#iterate through each column (bp position)
my $i=1;
while($i<$seqlen)
	{
	  $ncount=0;
	  #iterate through each entry
	  foreach my $seq ($aln->each_seq)
		{
		  my $res=$seq->subseq($i,$i);
		  if ($res eq "N"){$ncount++}
		}
	  if (($ncount/$numsamples) ge ($thresh/100))
		{
		  #remove column and decrease counter
		  my @tmp=($i-1,$i-1);
		  my $poss=\@tmp;
		  $aln=$aln->remove_columns($poss);
		  $seqlen=$aln->length;
		}
		else {$i++}
	}


#basic stats
$numsamples=$aln->num_sequences;
print "number of samples is now $numsamples\n";
$seqlen=$aln->length;
print "sequence length is now $seqlen\n";

#print to file
my $outfile=$infasta . $thresh;
print "filtered fasta writing to $outfile\n";
my $out = Bio::AlignIO->new(-file => ">$outfile", -format => 'fasta');
        $out->write_aln($aln);




























