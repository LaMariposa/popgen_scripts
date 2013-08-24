#modified to return numerator and denominator of Fst
#modified by Megan Supple, 2012
#
# BioPerl module for Bio::PopGen::PopStats
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code



# Let the code begin...


#package Bio::PopGen::PopStatsModified;
package PopStatsModified;
use strict;

use base qw(Bio::Root::Root);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($haploid) = $self->_rearrange([qw(HAPLOID)],@args);
  if( $haploid ) { $self->haploid_status(1) }
  return $self;
}



sub haploid_status{
    my $self = shift;
    return $self->{'haploid_status'} = shift if @_;
    return $self->{'haploid_status'};
}


# Implementation provided my Matthew Hahn, massaged by Jason Stajich


#' make emacs happy here
sub Fst {
   my ($self,$populations,$markernames) = @_;

   if( ! defined $populations || 
       ref($populations) !~ /ARRAY/i ) { 
       $self->warn("Must provide a valid arrayref for populations");
       return;
   } elsif( ! defined $markernames ||
	    ref($markernames) !~ /ARRAY/i ) {
       $self->warn("Must provide a valid arrayref for marker names");
       return;
   }
   my $num_sub_pops          = scalar @$populations;

   if( $num_sub_pops < 2 ) {
       $self->warn("Must provide at least 2 populations for this test, you provided $num_sub_pops");
       return;
   }

   # This code assumes that pop 1 contains at least one of all the
   # alleles - need to do some more work to insure that the complete 
   # set of alleles is seen.
   my $Fst;
   my ($TS_sub1,$TS_sub2);

   foreach my $marker ( @$markernames ) {
       # Get all the alleles from all the genotypes in all subpopulations
       my %allAlleles;
       foreach my $allele ( map { $_->get_Alleles() } 
			    map { $_->get_Genotypes($marker) } @$populations ){
	   $allAlleles{$allele}++;
       }
       my @alleles = keys %allAlleles;

       foreach my $allele_name ( @alleles ) {
	   my $avg_samp_size         = 0; # n-bar
	   my $avg_allele_freq       = 0; # p-tilda-A-dot

	   my $total_samples_squared = 0; # 
	   my $sum_heterozygote      = 0;

	   my @marker_freqs;

	   # Walk through each population, get the calculated allele frequencies
	   # for the marker, do some bookkeeping


	   foreach my $pop ( @$populations ) {
	       my $s = $pop->get_number_individuals($marker);

	       $avg_samp_size += $s;
	       $total_samples_squared += $s**2;

	       my $markerobj = $pop->get_Marker($marker);
	       if( ! defined $markerobj ) { 
		   $self->warn("Could not derive Marker for $marker ".
			       "from population ". $pop->name);
		   return;
	       }

	       my $freq_homozygotes = 
		   $pop->get_Frequency_Homozygotes($marker,$allele_name);
	       my %af = $markerobj->get_Allele_Frequencies();
	       my $all_freq = ( ($af{$allele_name} || 0));

	       $avg_allele_freq += $s * $all_freq;
	       $sum_heterozygote += (2 * $s)*( $all_freq - $freq_homozygotes);

	       push @marker_freqs, \%af;
	   }
	   my $total_samples =  $avg_samp_size;	# sum of n over i sub-populations
	   $avg_samp_size /= $num_sub_pops;
	   $avg_allele_freq /= $total_samples;
	   
	   # n-sub-c
	   my $adj_samp_size = ( 1/ ($num_sub_pops - 1)) *
	       ( $total_samples - ( $total_samples_squared/$total_samples));

	   my $variance              = 0; # s-squared-sub-A
	   my $sum_variance          = 0;
	   my $i = 0;		# we have cached the marker info
	   foreach my $pop ( @$populations ) {
	       my $s = $pop->get_number_individuals($marker);
	       my %af = %{$marker_freqs[$i++]};
	       $sum_variance += $s * (( ($af{$allele_name} || 0) - 
					$avg_allele_freq)**2);
	   }
	   $variance = ( 1 / (( $num_sub_pops-1)*$avg_samp_size))*$sum_variance;

	   # H-tilda-A-dot
	   my $freq_heterozygote = ($sum_heterozygote / $total_samples);

	   if( $self->haploid_status ) {
	       # Haploid calculations

	       my $T_sub1 = $variance - 
		   ( ( 1/($avg_samp_size-1))*
		     ( ($avg_allele_freq*(1-$avg_allele_freq))-
		       ( (($num_sub_pops-1)/$num_sub_pops)*$variance)));
	       my $T_sub2 = ( (($adj_samp_size-1)/($avg_samp_size-1))*
			      $avg_allele_freq*(1-$avg_allele_freq) ) +
			      ( 1 + ( (($num_sub_pops-1)*
				       ($avg_samp_size-$adj_samp_size))/ 
				      ($avg_samp_size - 1))) * 
				      ($variance/$num_sub_pops);


	       #to get total Fst from all alleles (if more than two) or all
	       #loci (if more than one), we need to calculate $T_sub1 and
	       #$T_sub2 for all alleles for all loci, sum, and then divide
	       #again to get Fst.
	       $TS_sub1 += $T_sub1;
	       $TS_sub2 += $T_sub2;

	   } else {
	       my $S_sub1 = $variance - ( (1/($avg_samp_size-1))*
					  ( ($avg_allele_freq*
					     (1-$avg_allele_freq)) - 
					    ((($num_sub_pops-1)/$num_sub_pops)*
					     $variance)-0.25*$freq_heterozygote ) );
	       my $S_sub2 = ($avg_allele_freq*(1-$avg_allele_freq)) - 
		   ( ($avg_samp_size/($num_sub_pops*($avg_samp_size-1)))*
		     ( ((($num_sub_pops*($avg_samp_size- $adj_samp_size))/
			 $avg_samp_size)*$avg_allele_freq*
			(1-$avg_allele_freq)) - 
		       ( (1/$avg_samp_size)* (($avg_samp_size-1)+
					      ($num_sub_pops-1)*
					      ($avg_samp_size-
					       $adj_samp_size) )*$variance ) - 
		       ( (($num_sub_pops*($avg_samp_size-$adj_samp_size))/
			  (4*$avg_samp_size*$adj_samp_size))*
			 $freq_heterozygote ) ) );

	       my $S_sub3 = ($adj_samp_size/(2*$avg_samp_size))*
		   $freq_heterozygote;

	       #Again, to get the average over many alleles or many loci,
	       #we will have to run the above for each and then sum the $S
	       #variables and recalculate the F statistics 
	       $TS_sub1 += $S_sub1;
	       $TS_sub2 += $S_sub2;
	   } 
       }
   }
   # $Fst_diploid = $S_sub1/$S_sub2;
   #my $Fit_diploid = 1 - ($S_sub3/$S_sub2);
   #my $Fis_diploid = ($Fit_diploid-$Fst_diploid)/(1-$Fst_diploid);
   $Fst = $TS_sub1 / $TS_sub2;

   return ($TS_sub1, $TS_sub2, $Fst);
}

1;
