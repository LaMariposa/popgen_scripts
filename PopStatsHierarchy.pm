#modified to account for population substructure following Weir GDAII p184
#modified by Megan Supple 18 May 2012
#
# BioPerl module for Bio::PopGen::PopStats
# Copyright Jason Stajich


#package Bio::PopGen::PopStatsHierarchy;
package PopStatsHierarchy;
use strict;
use Data::Dumper;
use List::MoreUtils qw/ uniq /;

use base qw(Bio::Root::Root);

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($haploid) = $self->_rearrange([qw(HAPLOID)],@args);
  if( $haploid ) { $self->haploid_status(1) }
  return $self;
}


sub Fst 
	{
   	  my ($self,$populations,$markernames) = @_;

	  #test to see if inputs ok
   	  if( ! defined $populations || ref($populations) !~ /ARRAY/i ) 
		{ 
		  $self->warn("Must provide a valid arrayref for populations");
       		  return;
   		} 
		elsif( ! defined $markernames || ref($markernames) !~ /ARRAY/i ) 
		{
       		  $self->warn("Must provide a valid arrayref for marker names");
       		  return;
   		}
	  if( scalar @$populations < 2 ) 
		{
       		  $self->warn("Must provide at least 2 populations for this test");
       		  return;
   		}

		  
	  #get information on populations(location) and subpopulations(pheno)
		  #loop through each race (pop_pheno)
	  my @phenos_all; my @locations_all;
	  my %pop_info;
	  foreach my $pop ( @$populations )
		{
		  my $name=$pop->name;
		  my $loc=$pop->source;
		  my $pheno=$pop->description;

		  push (@phenos_all, $pheno);
		  push (@locations_all, $loc);

		  $pop_info{$name}=$pop;
		}
	  my @phenos=uniq @phenos_all;
	  my @locations=uniq @locations_all;


	  # ???This code assumes that pop 1 contains at least one of all the alleles - 
		#need to do some more work to insure that the complete set of alleles is seen.
	  my $theta_hat_subS=0; my $R2=0; my $R3=0;

	  #look at each marker 
   	  foreach my $marker ( @$markernames ) 
		{
		  #initialize a subpop size array (assumes 2 phenos)
		  my $n_pops=scalar @locations;	#r=number of populations
		  my $n_phenos=scalar @phenos;
		  my @n; my @n_i_dot; my $n_dot_dot=0; my $s_dot=0;

		  #get population size information, iterate through all races
			#assumes 2 subpopulations for all populations, should not effect calculations because just adds zeros
                  for (my $i=0;$i<$n_pops;$i++)
			{
                          for (my $j=0;$j<$n_phenos;$j++)
				{
				  my $pop_p=$pop_info{$locations[$i] . "_" . $phenos[$j]};
                                  my $pop_size=$pop_p->get_number_individuals($marker);
                                  $n[$i][$j] = $pop_size;
                                  $n_i_dot[$i]+=$pop_size;
				  $n_dot_dot+=$pop_size;
				  if ($pop_size>0){$s_dot++;}
				}
			}

		  #calculate nc1, nc2, nc3
		  my $nc1=0; my $nc2=0; my $nc3=0;
		  my $nc1_part=0; my $nc2_part=0; my $nc3_part=0;
                  for (my $i=0;$i<$n_pops;$i++)
                        {
                          for (my $j=0;$j<$n_phenos;$j++)
                                {
				  $nc1_part+=($n_dot_dot-$n_i_dot[$i])*$n[$i][$j]*$n[$i][$j]/($n_i_dot[$i]*$n_dot_dot);
				  $nc3_part+=$n[$i][$j]*$n[$i][$j]/$n_i_dot[$i];
				}
			  $nc2_part+=$n_i_dot[$i]*$n_i_dot[$i];
			}                         
		  $nc1=$nc1_part/($n_pops-1);
		  $nc2=($n_dot_dot-($nc2_part/$n_dot_dot))/($n_pops-1);
		  $nc3=($n_dot_dot-$nc3_part)/($s_dot-$n_pops);

       		  # Get all the alleles from all the genotypes in all subpopulations
	       	  my %allAlleles;
       		  foreach my $allele ( map { $_->get_Alleles() } map { $_->get_Genotypes($marker) } @$populations )
			{$allAlleles{$allele}++;}
 	     	  	my @alleles = keys %allAlleles;

		  #iterate through each allele
       		  foreach my $allele_name ( @alleles ) 
			{

			  # Walk through each population, get the calculated allele frequencies, etc
			  #loop through each population (i) and look at both subpopulations (j)
			  my @p_tilda_A_i_j; my @p_tilda_A_i_dot; my $p_tilda_A_dot_dot=0;
			  my @H_tilda_A_i_j;
			  for (my $i=0;$i<$n_pops;$i++)
				{
				  for (my $j=0;$j<$n_phenos;$j++)
					{
					  my $pop=$pop_info{$locations[$i] . "_" . $phenos[$j]};
				       	  my $markerobj = $pop->get_Marker($marker);
			 	       	  if( ! defined $markerobj ) {$self->warn("Could not derive Marker for $marker "."from population ". $pop->name);return;}

					  $H_tilda_A_i_j[$i][$j]=$pop->get_Frequency_Heterozygotes($marker,$allele_name);

			       	  	  my %af = $markerobj->get_Allele_Frequencies();
	       				  $p_tilda_A_i_j[$i][$j] = ( ($af{$allele_name} || 0));
	       				  $p_tilda_A_i_dot[$i]+=$n[$i][$j]*$p_tilda_A_i_j[$i][$j];
					  $p_tilda_A_dot_dot+=$n[$i][$j]*$p_tilda_A_i_j[$i][$j];
					}
				  $p_tilda_A_i_dot[$i]/=$n_i_dot[$i];
	   			}
			  $p_tilda_A_dot_dot/=$n_dot_dot;
	
			  #calc MSP, MSS, MSI, and MSG
			  my $MSP=0; my $MSS=0; my $MSI=0; my $MSG=0;
	                  for (my $i=0;$i<$n_pops;$i++)
        	                {
                	          for (my $j=0;$j<$n_phenos;$j++)
                        	        {
					  $MSS+=$n[$i][$j]*(($p_tilda_A_i_j[$i][$j]-$p_tilda_A_i_dot[$i])**2);
					  $MSI+=$n[$i][$j]*($p_tilda_A_i_j[$i][$j]*(1-$p_tilda_A_i_j[$i][$j])-($H_tilda_A_i_j[$i][$j]/4));
					  $MSG+=$n[$i][$j]*$H_tilda_A_i_j[$i][$j];
					}
				  $MSP+=$n_i_dot[$i]*($p_tilda_A_i_dot[$i]-$p_tilda_A_dot_dot);
				}
			  $MSP*=2/($n_pops-1);
			  $MSS*=2/($s_dot-$n_pops);
			  $MSI*=2/($n_dot_dot-$s_dot);
			  $MSG*=1/(2*$n_dot_dot);
			  #calculate R3, R2, and theta_hat_s for the alleles
			  my $R3_temp=(($MSP-$MSI)/(2*$nc2))+((($nc2-$nc1)*($MSS-$MSI))/(2*$nc2*$nc3));
			  my $R2_temp=(($MSP-$MSI)/(2*$nc2))+((($nc2-$nc1)*($MSS-$MSI))/(2*$nc2*$nc3))+(($MSI+$MSG)/2);
			  #average R3, R2 over alleles and loci by summing
			  $R3+=$R3_temp;
			  $R2+=$R2_temp;
			}	
   		}
	  $theta_hat_subS=$R3/$R2;
   	  return ($R3, $R2, $theta_hat_subS);
	}

1;
