popgen_scripts
==============
Scripts for various population genomics analyses:

fasta2popgen.pl
---------------
Calculates various population genetic metrics from aligned fasta entries. Measures include Fst, hierarchical Fst, genotype by phenotype association, heterozygosity, pi, segregating sites, Tajima's D, D*, singletons, and SNP type (private vs. shared). Input is an aligned fasta file and a text file with sample information (population, phenotype, sample identifier). Output are population genomic metrics per position. Requires BioPerl, FormatGenos.pm (from the genome resequencing repository), and other modules from this repository (PopGen.pm, PopStatsHierarchy.pm, PopStatsModified.pm).

fst_stats.pl
------------
Calculates sliding window Fst and baseline Fst with jackknife confidence intervals. Input is per position Fst formatted as output from fasta2popgen.pl script.

sliding_av.pl
-------------
Calculates sliding window average.  Input file has columns: contig, position, number of samples with data for population 1, number of samples with data for population 2, value.  Input parameters are window size, step size, required minimum proportion of sites with data in each window, required minimum proportion of individuals of each population for each site.  Output is a text file with sliding window values.
 
