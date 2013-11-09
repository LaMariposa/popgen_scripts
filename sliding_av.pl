#!/usr/bin/perl

#sliding_av.pl by Megan Supple
#created 24 Aug 2013
#last modified 24 Aug 2013
#
#script to calculate sliding window average
#usage sliding_av.pl InputFile.txt WindowSize StepSize MinPropSites MinPropIndivs
#	InputFile.txt is formatted like the dxy.txt output file from fasta2popgen.pl
#		contig, position, pop1#genoed, pop2#genoed, dxy
#	WindowSize is the size of the sliding window
#	StepSize is the size of the step when moving the sliding window
#	MinPropSites is the minimum proportion of markers needed to calculate average at a window
#	MinPropIndivs is the minimum proportion of individuals of each population needed to use a site



use lib $ENV{PERL5LIB};
 
use strict; 
use warnings;
use Getopt::Long;
use Data::Dumper;

my $usage = "Usage: sliding_av.pl <InputFile.txt> <WindowSize> <StepSize> <MinPropSites> <MinPropIndivs>
arguments (required):
	<InputFile.txt> is per position dxy from fasta2popgen.pl
	<WindowSize> is the size of the sliding window
	<StepSize> is the size of the step when moving the sliding window
	<MinPropSites> is the minimum proportion of markers needed to calculate average at a window
	<MinPropIndivs> is the minimum proportion of individuals of each population needed to use a site
";

die "$usage" unless (@ARGV == 5);

#read in command line arguments
my ($infile,$windowSize,$stepSize,$minProp, $minPropInds)=@ARGV;

#open input file and output file for sliding window
open(INFILE, $infile)||die "can't open input file. $!\n";
my $outfile="sliding_av_" . $windowSize . "x" . $stepSize . "x" . $minProp . "x" . $minPropInds .".txt";
open(OUTFILE, ">$outfile");

#print input parameters to output file
print OUTFILE "####################\n";
print OUTFILE "#input file=$infile\n";
print OUTFILE "#window size=$windowSize\n";
print OUTFILE "#step size=$stepSize\n";
print OUTFILE "#minimum proportion of sites=$minProp\n";
print OUTFILE "#minimum proportion of individual=$minPropInds\n";

#read header info from input file and write to output file
my $pop1n; my $pop2n;

my $line=<INFILE>; print OUTFILE "$line";#read in header border
$line=<INFILE>; print OUTFILE "$line";	#read in population 1 sample size
	chomp $line;
	my @entry=split("=",$line);
	$pop1n=$entry[1];
$line=<INFILE>; print OUTFILE "$line";	#read in population 1 sample list
$line=<INFILE>; print OUTFILE "$line";	#read in population 2 sample size
        chomp $line;
	@entry=split("=",$line);
        $pop2n=$entry[1];
$line=<INFILE>; print OUTFILE "$line";	#read in population 2 sample list
$line=<INFILE>; print OUTFILE "$line";	#read in header border
$line=<INFILE>;				#read column header

print "n1=$pop1n and n2=$pop2n\n";
print OUTFILE "contig\tposition\taverage\n";

print "calculating sliding window\n";

#declare variables to track window
my $window_start=1;
my $window_end=$windowSize;

#declare other variable
my $contig;
my $av;
my $mid_pos;
my @current_set=();	#an array of arrays holding the current set of markers

#process input file until EOF
while($line=<INFILE>)
	{
	  #break up line into component parts
	  my @entry=split(" ", $line);
	  $contig=$entry[0];
	  my $pos=$entry[1];

	  #calculate sliding window average
	  #determine if current entry is in the window
	  if ($pos>=$window_start && $pos<=$window_end)
		{
		  #marker is in window
		  #if enough individuals, add it to current marker set
		  if($entry[2]>=$minPropInds*$pop1n && $entry[3]>=$minPropInds*$pop2n)
			{
			  push (@current_set, [@entry]);
			}
		}
	  else
		{
	  	  while ($pos<$window_start || $pos>$window_end)
			{
			  #marker is not in window, need to calculate average windows until it is in window		  			
		  	  #calculate average from previous window if there are enough markers, otherwise print NAs
			  if (@current_set<$windowSize*$minProp)
				{
				  $mid_pos=int($window_start+.5*($windowSize-1));
                        	  print OUTFILE "$contig\t$mid_pos\tNA\n";
				}
			  else
				{
				  #calc average
				  my $sum=0;
				  my $n=0;
				  for (my $i=0;$i<@current_set;$i++)
					{
					  $sum+=$current_set[$i][4];
					  $n++;
					}
				  $av=eval{$sum/$n};

				  #calc mid position
				  $mid_pos=int($window_start+.5*($windowSize-1));
				  #print results to outfile	
				  if (defined $av) {print OUTFILE "$contig\t$mid_pos\t$av\n";}
					else {print OUTFILE "$contig\t$mid_pos\tNA\n";}
				}	
 
			  #reset current window
			  $window_start+=$stepSize;
			  $window_end=$window_start+$windowSize-1;
			  #reset current marker set by removing elements that are not in new window
			  for (my $i=0;$i<@current_set;$i++)
				{
				  #look at each element to see if it is in new window
				  if ($current_set[$i][1]<$window_start || $current_set[$i][1]>$window_end)
					{
					  #not in new window so remove from current set
					  splice(@current_set, $i, 1); 					  
					  #removed an entry so need to reindex
					  $i--; 
					}
				}
			}

	  	  #marker is now in current window so push into current array if enough samples
                  if($entry[2]>=$minPropInds*$pop1n && $entry[3]>=$minPropInds*$pop2n)
                        {
                          push (@current_set, [@entry]);
                        }
		}
	}


#after EOF calculate final window if there are enough markers
if (@current_set<$windowSize*$minProp)
	{
          $mid_pos=int($window_start+.5*($windowSize-1));
          print OUTFILE "$contig\t$mid_pos\tNA\n";
	}
else
        {
          #calculate average
          my $sum=0;
          my $n=0;
          for (my $i=0;$i<@current_set;$i++)
  		{
                  $sum+=$current_set[$i][4];
		  $n++;
                }
          $av=eval{$sum/$n};

          #calc mid position
          $mid_pos=int($window_start+.5*($windowSize-1));
          #print results to outfile
          if (defined $av) {print OUTFILE "$contig\t$mid_pos\t$av\n";}
		else {print OUTFILE "$contig\t$mid_pos\tNA\n";}
      	}



close INFILE;
close OUTFILE;
print "done!\n";










 












