#perl modules for calculating population genetic metrics
#Megan Supple
#created Nov 2011
#last modified 14 Aug 2013

package PopGen;

use strict;
use warnings;

use Bio::PopGen::Individual;
use Bio::PopGen::Statistics;
use List::Util qw(max sum);
use Text::NSP::Measures::2D::Fisher::twotailed;
use Data::Dumper;

#####################################################################################
#Fst
#calculates Fst at each position along genomic intervals from genotype objects
#input is pointers to population names, popoulation data, contig name, contig size, header information
#ouputs text file with Fst per position--contig, position, # genotyped in pop 1, #genotyped in pop1, s1, s2, Fst

sub Fst
{ 
	use PopStatsModified;
	
	#read in inputs
	my ($pop_names_p, $pops_p, $contig_p, $size_p, $header_p)=@_;

	#open output file
	open(OUT,">fst.txt"); #fst output file
	#print header information
	print OUT @$header_p;
	print OUT "contig\tposition\tn_pop1\tn_pop2\ts1\ts2\tFst\n";
	
	#calculate fst
	#for each position calculate fst and print result to file
	my $pos=1;
	while ($pos<=$$size_p)
		  	{
	  		  #determine number of individual genotyped 
			  my @genoed; #array of sample sizes indexed by population
			  #check each population
			  for (my $j=0;$j<scalar @$pop_names_p;$j++)
				{
				  #get an array of individuals that are genotyped at the marker	  
				  my @genoed_inds = \$$pops_p[$j]->get_Individuals(-marker => $pos);
				  $genoed[$j]=scalar @genoed_inds;
		 		}

			  #calculate fst if any genotypes
			  if ($genoed[0]>0 && $genoed[1]>0)
				{
				  my $stats = PopStatsModified->new();
				  #maker marker into an array to make fst calculator happy
				  my @marker;
				  $marker[0]=$pos;
				  #use eval to allow recovery from divide by zero and skip window
				  my ($s1, $s2, $fst)=eval {$stats->Fst(\@$pops_p,\@marker)};
				  #if able to calculate fst, print result to file
				  if (defined $fst) {print OUT "$$contig_p\t$pos\t$genoed[0]\t$genoed[1]\t$s1\t$s2\t$fst\n"}
					else {print OUT "$$contig_p\t$pos\t$genoed[0]\t$genoed[1]\t0\t0\tNA\n"}
				}	  
		  	  #move to start of next position
 		  	  $pos++;
	   	 	}
 
 	close OUT;
}



#####################################################################################
#association
#calculates genotype by phenotype association for each position in the contig
#input is pointers to population names, popoulation data, contig name, contig size, header information
#ouputs text file with association per position--contig, position, # genotyped in pop 1, #genotyped in pop1,  pval of association, if there perfect association

sub association
{ 
	#read in inputs
	my ($pop_names_p, $pops_p, $contig_p, $size_p, $header_p)=@_;

	#open output file
	open(OUT2,">assoc.txt");
	#print header information 
	print OUT2 @$header_p;
	print OUT2 "contig\tposition\tn_pop1\tn_pop2\tpval\tis_fixed\n";
	
	#calculate genotype by phenotype association
	#for each position calculate and print result to file
	my $pos=1;
	while ($pos<$$size_p)
	  	{	
		  #determine number of individuals genotyped
		  my @genoed; ##array of sample sizes indexed by population
		  #check each population	
	  	  for (my $j=0;$j<@$pop_names_p;$j++)  #j tracks the population
			{
			  #get an array of individuals that are genotyped at the marker	  
			  my @genoed_inds = \$$pops_p[$j]->get_Individuals(-marker => $pos);
			  $genoed[$j]=scalar @genoed_inds;
			}

	 	  #calc genotype by phenotype association
		  #initialize data stucture
		  my @allele_counts;
		  for (my $a=0;$a<@$pop_names_p;$a++)
			{
			  $allele_counts[$a]{'A'}=0;
			  $allele_counts[$a]{'T'}=0;
			  $allele_counts[$a]{'C'}=0;
			  $allele_counts[$a]{'G'}=0;
			}
		  #for each population, for each individual get phenotype and genotype
		  for (my $j=0;$j<@$pop_names_p;$j++)  #j tracks the population
	 		{
			  #get genotypes at the position for each individual in the population 
			  my @genotypes=\$$pops_p[$j]->get_Genotypes(-marker => $pos);
			  #for each individual in the population
	 		  my @genoed_inds = \$$pops_p[$j]->get_Individuals(-marker => $pos);
	 		  for (my $k=0;$k<scalar @genoed_inds;$k++) #k tracks the individual
				{
				  #for each allele
				  for (my $l=0;$l<2;$l++) #l tracks the allele
					{
					  #get genotype and add to allele count
					  $allele_counts[$j]{"${$${$genotypes[$k]}{_alleles}}[$l]"}++;
					}
				}
			}

 		  #calculate allele counts in super population
		  #get allele counts for all populations
		  #initalize hash
		  my %total_af = ( 'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0 );
		  for (my $j=0;$j<@$pop_names_p;$j++) 
			{
			  $total_af{'A'}+=$allele_counts[$j]{'A'};  
			  $total_af{'T'}+=$allele_counts[$j]{'T'};
			  $total_af{'C'}+=$allele_counts[$j]{'C'};
			  $total_af{'G'}+=$allele_counts[$j]{'G'};
			}

		  #determine if triallelic and shows variation
		  my $n_alleles=0;
      	          if($total_af{'A'}>0) {$n_alleles+=1}
           	  if($total_af{'T'}>0) {$n_alleles+=1}
               	  if($total_af{'C'}>0) {$n_alleles+=1}
       		  if($total_af{'G'}>0) {$n_alleles+=1}

  	 	  #calc assoc if not triallelic and does show variation		  
		  if($n_alleles==2)
			{
			  #determine major allele
			  my $max = max values %total_af;
			  my $major= +{ reverse %total_af }->{$max};
			
			  #calculate total number of alleles
			  my $n_alleles=sum(values %total_af);

			  #determine if there is variation--already did above
			  #if(($max/$n_alleles)!=1) 
			 	#{
				  #there is variation so calc fishers
				  #now assume just 2 populations
				  #total number of alleles
				  my $npp=$n_alleles;  
				  #total number of major allele
				  my $np1=$total_af{$major};
				  #total number of alleles in pop1
				  my $n1p=$allele_counts[0]{'A'}+$allele_counts[0]{'T'}+$allele_counts[0]{'C'}+$allele_counts[0]{'G'};
				  #number of major allele in pop1
				  my $n11=$allele_counts[0]{$major};
				  my $twotailed = calculateStatistic( n11=>$n11,
                               		n1p=>$n1p,
            	              		np1=>$np1,
      	                      		npp=>$npp);
		  		  #test if fixed difference
				  my $fixed_flag;
				  if ($n11==0 && $npp-$n1p-$np1==0){$fixed_flag="fixed"}
					elsif ($n1p==$n11 && $np1==$n11) {$fixed_flag="fixed"}
					else {$fixed_flag="no"}
				  
				  #print out results
				  print OUT2 "$$contig_p\t$pos\t$genoed[0]\t$genoed[1]\t$twotailed\t$fixed_flag\n";
				#}		  	
	  		}
	 	  	  #move to next position
  		  	  $pos++;	  			
		}
	  close OUT2;
}		

	
#####################################################################################
#Fst_3level
#caluculate 3 level hierarchy based on Weir GDAII
#calculates Fst at each position
#input is pointers to population names, popoulation data, contig name, contig size, header information
#ouputs text file with Fst along intervals--contig, position, R2, R3, theta_hat_sub_s

sub Fst_3level
{
use PopStatsHierarchy;
        #read in inputs
        my ($pop_names_p, $pops_p, $contig_p, $size_p, $header_p)=@_;

        #open output file
        open(OUT,">fst3level.txt"); #fst output file
	#print header information
        print OUT @$header_p;

	#calculate the population size of each phenotype
	my %pheno_size;
        my $n_pops=scalar @$pop_names_p;
        for (my $q=0;$q<$n_pops;$q++)
                {
                  my $num_inds = \$$pops_p[$q]->get_number_individuals;
		  #print OUT "$$num_inds individuals in $$pop_names_p[$q]\n";
		  #determine phenotype
		  my $pheno=\$$pops_p[$q]->description;
		  #my @info=split(/_/,$$pop_names_p[$q]);
		  if ($pheno_size{$$pheno}){$pheno_size{$$pheno}+=$$num_inds;}
			else{$pheno_size{$$pheno}=$$num_inds;}
                }
	my @phenos=keys %pheno_size; 
#	while (my ($pheno,$size)=each(%pheno_size)) {print OUT "number of individuals of $pheno type is $size\n";}
        print OUT "contig\tposition\tn_$phenos[0]\tn_$phenos[1]\tR3\tR2\ttheta_hat_subS\n";

        #calculate fst 3 level heirarchy for each position and print result to file
        my $pos=1;
        while ($pos<=$$size_p)
		{
                  #determine number of individuals genotyped
		  my @genoed; #array of sample sizes indexed by phenotype
                  #check each population
                  for (my $j=0;$j<scalar @$pop_names_p;$j++)
                	{
                          #get an array of individuals that are genotyped at the marker
                          my @genoed_inds = \$$pops_p[$j]->get_Individuals(-marker => $pos);
			  my $n=scalar @genoed_inds;
			  #get phenotype
			  my $pheno= \$$pops_p[$j]->description;
			  if ($$pheno eq $phenos[0]){$genoed[0]+=$n;}
			    elsif ($$pheno eq $phenos[1]){$genoed[1]+=$n;}
			    else {print "uh oh, new phenotype\n";}				
                        }
		  #calculate fst if any individuals in each phenotype
                  if ($genoed[0]>0 && $genoed[1]>0)
			{
                          #calc fst
                          my $stats = PopStatsHierarchy->new();
                          #maker marker into an array to make fst calculator happy
                          my @marker;
                          $marker[0]=$pos;
                          #use eval to allow recovery from divide by zero and skip window
                          my ($r3, $r2, $theta_hat_subS)=eval {$stats->Fst(\@$pops_p,\@marker)};
                          #if able to calculate fst, print result to file
                          if (defined $theta_hat_subS) {print OUT "$$contig_p\t$pos\t$genoed[0]\t$genoed[1]\t$r3\t$r2\t$theta_hat_subS\n"}
                            else {print OUT "$$contig_p\t$pos\t$genoed[0]\t$genoed[1]\t0\t0\tNA\n"}
                        }
                   #move to start of next position
                   $pos++;
                 }
        close OUT;
}


############################################################################################
#selection
#calculate various measures used to test for selection for each population and all combined
#metrics: heterozygosity, pi, segregating sites, Tajima's D, D*, singletons, SNP type (private vs shared)
#input is pointers to population names, popoulation data, contig name, contig size, header information
#putputs a file with all metrics

sub selection
{
	my @alleles=([0,0,0,0],[0,0,0,0]);
        #read in inputs
        my ($pop_names_p, $pops_p, $contig_p, $size_p, $header_p)=@_;

        #open output file
        open(OUT,">selection.txt"); 
        #print header information
	print OUT @$header_p;
        print OUT "contig\tposition\t";
	for (my $q=0;$q<scalar @$pop_names_p;$q++)
                {
                  print OUT "$$pop_names_p[$q]_size\t$$pop_names_p[$q]_het\t$$pop_names_p[$q]_pi\t$$pop_names_p[$q]_segsites\t";
		  print OUT "$$pop_names_p[$q]_tajD\t$$pop_names_p[$q]_dstar\t$$pop_names_p[$q]_singleton\t";
                }
	print OUT "het_all\tpi_all\tsegsites_all\ttajD_all\tdstar_all\tsingletons_all\tsnp_type\n";

        #calculate selection metrics
	#for each position calculate metrics and print to file 
        my $pos=1;
        while ($pos<$$size_p)
                {
		  my @temp_inds_all; #array of individuals from all populations
		  my $het_all=0;	#count of heterozygotes in all samples
		  my @alleles=([0,0,0,0],[0,0,0,0]);
		  #print contig position info
		  print OUT "$$contig_p $pos\t";
                  #for each population, for each individual get phenotype and genotype
                  for (my $j=0;$j<scalar @$pop_names_p;$j++)  #j tracks the population
                	{
                       	  #get genotypes at the position for each individual in the population
                          my @genotypes=\$$pops_p[$j]->get_Genotypes(-marker => $pos);
                          #for each individual in the population
                          my @genoed_inds = \$$pops_p[$j]->get_Individuals(-marker => $pos);
			  #get the number of individuals genotyped in that population
			  my $n_pop=scalar @genoed_inds;
			  print OUT "$n_pop\t";

		  	  #rebuild info for just single marker if any individuals genoed
			  my @temp_inds;
			  my $het_pop=0;
		  	  if ($n_pop==0){print OUT "0\tNA\t0\tNA\tNA\t0\t";}
			 	  else{
                                        for (my $k=0;$k<scalar @genoed_inds;$k++) #k tracks the individual
                                                {
						  my @temp_genotype = ${$genoed_inds[$k]}->get_Genotypes(-marker => $pos);
						  my $ind = Bio::PopGen::Individual->new(-unique_id => $genoed_inds[$k], -genotypes => \@temp_genotype);
						  push @temp_inds, $ind;		
						  push @temp_inds_all, $ind;
						
						  #count heterozygotes
						  if (${${$temp_genotype[0]}{_alleles}}[0] ne ${${$temp_genotype[0]}{_alleles}}[1])
							{
							  $het_all++;
							  $het_pop++;
							}

						  #track alleles in the population
						  if (${${$temp_genotype[0]}{_alleles}}[0] eq "A"){$alleles[$j][0]=1;}
                                                  if (${${$temp_genotype[0]}{_alleles}}[0] eq "T"){$alleles[$j][1]=1;}
                                                  if (${${$temp_genotype[0]}{_alleles}}[0] eq "C"){$alleles[$j][2]=1;}
                                                  if (${${$temp_genotype[0]}{_alleles}}[0] eq "G"){$alleles[$j][3]=1;}
                                                  if (${${$temp_genotype[0]}{_alleles}}[1] eq "A"){$alleles[$j][0]=1;}
                                                  if (${${$temp_genotype[0]}{_alleles}}[1] eq "T"){$alleles[$j][1]=1;}
                                                  if (${${$temp_genotype[0]}{_alleles}}[1] eq "C"){$alleles[$j][2]=1;}
                                                  if (${${$temp_genotype[0]}{_alleles}}[1] eq "G"){$alleles[$j][3]=1;}
					        }
					      
					#print hets
					print OUT "$het_pop\t";
					my $stats = Bio::PopGen::Statistics->new();
					#calculate pi
                                      	my $pi = $stats->pi(\@temp_inds);
                                        print OUT "$pi\t";
					#calculate segregating sites
					my $segsites = $stats->segregating_sites_count(\@temp_inds);
					print OUT "$segsites\t";
					#calculate tajima's D
					my $tajd = eval{$stats->tajima_D(\@temp_inds)};
					if ($tajd){print OUT "$tajd\t";}
						else {print OUT "NA\t";}
					#calculate Fu & Li's D* and F*
					my $d_star=eval{$stats->fu_and_li_D_star(\@temp_inds)};
                                        if ($d_star){print OUT "$d_star\t";}
                                            else {print OUT "NA\t";}
					#my $f_star;
                                        #if($n_pop>1){$f_star=eval{$stats->fu_and_li_F_star(\@temp_inds)};}
                                        #if ($f_star){print OUT "fstar=$f_star\t";}
                                        #        else {print OUT "fstar=NA\t";}
					#calculate singletons
					my $singletons=$stats->singleton_count(\@temp_inds);
					print OUT "$singletons\t";
				     }
                        }

		  #calculate diversity measures for whole population
		  if (scalar @temp_inds_all==0){print OUT "0\tNA\t0\tNA\tNA\t0\tuntyped\n";}
			  else
				{
				   print OUT "$het_all\t";
				   my $stats = Bio::PopGen::Statistics->new();
                                   my $pi = $stats->pi(\@temp_inds_all);
                        	       	   print OUT "$pi\t";
				   my $segsites = $stats->segregating_sites_count(\@temp_inds_all);
                                           print OUT "$segsites\t";
				   my $tajd = eval{$stats->tajima_D(\@temp_inds_all)};
				   	if ($tajd){print OUT "$tajd\t";}
                                           else {print OUT "NA\t";}
                                   my $d_star=eval{$stats->fu_and_li_D_star(\@temp_inds_all)};
                                        if ($d_star){print OUT "$d_star\t";}
                                           else {print OUT "NA\t";}
				   #my $f_star;
                                        #if(@temp_inds_all>1){$f_star=eval{$stats->fu_and_li_F_star(\@temp_inds_all)};}
                                        #if ($f_star){print OUT "fstar=$f_star\t";}
                                        #        else {print OUT "fstar=NA\t";}
				   my $singletons=$stats->singleton_count(\@temp_inds_all);
                                            print OUT "$singletons\t";

				   #determine snp type
				   my @allele_tally=(0,0,0,0); my $pop0_tally=0; my $pop1_tally=0;
				   for (my $i=0;$i<4;$i++)
					{
					  $allele_tally[$i]=$alleles[0][$i]+$alleles[1][$i];
					  $pop0_tally+=$alleles[0][$i];
					  $pop1_tally+=$alleles[1][$i];
					}
				   #count alleles seen 
				   my $num_alleles=0;
				   for (my $i=0;$i<4;$i++){if ($allele_tally[$i]>0){$num_alleles++}}
				   if ($num_alleles==1){print OUT "invariant\n";}
				   if ($num_alleles>2){print OUT "multiallelic\n";}
				   if ($num_alleles==2)
					{
					  if ($allele_tally[0]==1 || $allele_tally[1]==1 || $allele_tally[2]==1 ||$allele_tally[3]==1)
						{
						  if ($allele_tally[0]+$allele_tally[1]+$allele_tally[2]+$allele_tally[3]==2)
							{print OUT "fixed\n";}
				    		    else 
							{
							  #determine which population it is private to
							  if ($pop0_tally==2){print OUT "private_1\n";}
								  elsif ($pop1_tally==2){print OUT "private_2\n";}
								  else {print OUT "oops2\n";}
							}
					}
					  elsif ($allele_tally[0]+$allele_tally[1]+$allele_tally[2]+$allele_tally[3]==4)
						{print OUT "shared\n"}
					  else {print OUT "oops\n"}

				}
			}	
		  #move to next position	
                  $pos++;
               }
         close OUT2;
}



#####################################################################################
#dxy
#caluculate absolute genetic divergence (dxy) from Nei 1987, eq 10.20
#calculates dxy at each position
#input is pointers to population names, popoulation data, contig name, contig size, header information
#ouputs text file with dxy along intervals--contig, position, population 1 samples size, population 2 sample size, dxy

sub dxy
{
        
	#read in inputs
        my ($pop_names_p, $pops_p, $contig_p, $size_p, $header_p)=@_;

	#open output file
	open(DXY, ">dxy.txt");
	#print out header information
	print DXY @$header_p;
	print DXY "contig\tposition\tn_pop1\tn_pop2\tdxy\n";

	#calculate dxy for each position and print to an output file
	my $pos=1;
	while ($pos<$$size_p)
		{
		  #determine number of individuals genotyped
		  my @genoed; ##array of sample sizes indexed by population
		  #check each population
		  for (my $j=0;$j<@$pop_names_p;$j++) #j tracks the population
			{
			  #get an array of individuals that are genotyped at the marker
			  my @genoed_inds = \$$pops_p[$j]->get_Individuals(-marker => $pos);
			  $genoed[$j]=scalar @genoed_inds;
			}


		  #get allele frequencies
		  #initialize data stucture
		  my @allele_counts;
		  for (my $a=0;$a<@$pop_names_p;$a++)
			{
			  $allele_counts[$a]{'A'}=0;
			  $allele_counts[$a]{'T'}=0;
			  $allele_counts[$a]{'C'}=0;
			  $allele_counts[$a]{'G'}=0;
			}
		  #for each population, for each individual get phenotype and genotype
		  for (my $j=0;$j<@$pop_names_p;$j++) #j tracks the population
			{
			  #get genotypes at the position for each individual in the population
			  my @genotypes=\$$pops_p[$j]->get_Genotypes(-marker => $pos);
			  #for each individual in the population
			  my @genoed_inds = \$$pops_p[$j]->get_Individuals(-marker => $pos);
			  for (my $k=0;$k<scalar @genoed_inds;$k++) #k tracks the individual
				{
				  #for each allele
				  for (my $l=0;$l<2;$l++) #l tracks the allele
					{
					  #get genotype and add to allele count
					  $allele_counts[$j]{"${$${$genotypes[$k]}{_alleles}}[$l]"}++;
					}
				}
			}

  		  #calculate allele counts in super population
		  #get allele counts for all populations
		  #initalize hash
		  my %total_af = ( 'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0 );
		  for (my $j=0;$j<@$pop_names_p;$j++)
			{
			  $total_af{'A'}+=$allele_counts[$j]{'A'};
			  $total_af{'T'}+=$allele_counts[$j]{'T'};
			  $total_af{'C'}+=$allele_counts[$j]{'C'};
			  $total_af{'G'}+=$allele_counts[$j]{'G'};
			}
		  #determine if triallelic and shows variation
		  my $n_alleles=0;
       		  if($total_af{'A'}>0) {$n_alleles+=1}
            	  if($total_af{'T'}>0) {$n_alleles+=1}
                  if($total_af{'C'}>0) {$n_alleles+=1}
        	  if($total_af{'G'}>0) {$n_alleles+=1}
		  
		  #if a biallelic SNP
		  if($n_alleles==2)
			{
print "bialleleic snp $$contig_p at $pos\n";
			  #print Dumper(@allele_counts);
			}
		  #move to next position
		  $pos++;
		}	
	close DXY;
}

1;
