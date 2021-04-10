#!/usr/bin/perl
use strict;
use warnings;
use Math::Round;
use Cwd;
#################################################################################
# Written by Nick Kent, Jan 19 2018
# Fiddled with: Nick Kent, Jan 19 2018 - tested
# Modified by Nick Kent, Aug 19 2019 - change to use cwd
#################################################################################
# USAGE:- perl SAM2dIEL_sgr_full.plx
#
# This script will take a SORTED Bowtie 1 SINGLE read alignment .sam format file from a CPSA exp
# and will generate chromosome-specific .sgr files containing 3MA-smoothed read mid-
# point frequency values for read ends as a proxy for MNase cleavage sites. This allows you
# to perform "digital Indirect End Labelling" or dIEL. This process is both subtly and fundamentally
# different to what SAM2PartN_sgr_full.plx does, so make sure you use the right script!
#
# This script will NOT run on a weedy Windows Laptop, but requires a high memory
# Linux machine. Typically you will need RAM eqivalent to the 0.75 * size 
# of your .sam file. There is also a limit to the number of reads that can be processed,
# and the number of genomic bins. These limits are both = 2^31 - the max size of a perl array. This
# script was designed for model eukaryote (yeast, Drosophila, Arabidopsis) and prokaryotic
# genomes, but may not work for mammalian size genomes.
#
# You must do a Bowtie 1-version SINGLE READ alignment of the Read1.fastq file from your paired end
# CPSA raw data. This file contains the highest quality reads, and the end of each read is equivalent to
# an MNase cleavage in the original DNA sample. 
# A typical bowtie CL might be:
# ./bowtie -n 0  -k 1 --best --trim3 39 -p 12 --sam indexes/S_cer_full_refseq WTyeastS5_R1_001.fastq WT_single36bp.sam
#
# User must then SORT the .sam file with SAMtools, e.g.:
# ./samtools sort Input.sam -@ 12 -O 'sam' -T Tempfilename -o Input_sorted.sam
#
# User MUST specify correctly the effective read length e.g. 36bp as above.
#
# SUGGESTION: once finished the user could "cat" (bash shell) the individual .sgr 
# files into one "whole genome" file, and use namechanger.plx to alter Chromosome IDs.
#
#
# NOTE: this script is the bastard love-child of Histo_Yeast_dIEL.plx and SAM2PartN_sgr_full.plx
#
# To do: ERROR TRAPPING (e.g. non-sorted .sam; general fuck-wittery;); rename test variables.
################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $indir_path   - The directory containing the SORTED .sam file to be processed
# $outdir_path  - The directory to store the .sgr output files
# $bin_width    - The histogram bin width value in bp (10bp is good for most things)
# $read_length 	- The length of the trimmed reads used in the Bowtie alignment in bp
################################################################################

my $inA_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my $bin_width = 20;
my $read_length = 36;

################################################################################
################################################################################
# MAIN PROGRAM - you should not need to edit below this line
################################################################################
################################################################################
print "Started:\t", `date`."\n";

# define some variables

my $infile_A;
my @line_A;
my @line_B;
my @files_A;
my @SAM_LN;
my @SAM_SN;
my @SAM_chr_id;
my @SAM_chr_length;
my @test_ID;
my @test_pos;
my $read_counter = 0;
my %SAM_map; 
my $SAM_data_count;
my $mapsize;
my $arguement_string;


################################################################################
# Find Chr IDs and gather specific data from .sam
################################################################################

# store input file name in an array
opendir(DIR,$inA_indir_path) || die "Unable to access file at: $inA_indir_path $!\n";

@files_A = readdir(DIR);

# process the input file within indir_path
foreach $infile_A (@files_A){    

    # ignore hidden files and only get those with the correct ending
    if (($infile_A !~ /^\.+/) && ($infile_A =~ /.*\.sam/)){
    
    
# While we're at it, let's print some useful information
print "Frequency distributions will be binned in $bin_width bp intervals \n";
print "\nRead length is set to:$read_length bp\n ";
print "\nFound, and processing, $infile_A \n\n";
print "Going to read SAM data into memory - this might take a while. \n";


open(IN, "$inA_indir_path/$infile_A")
            || die "Unable to open $infile_A: $!";
        
       
	
	# loop through top of infile to get header values
        while(<IN>){
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_A = split('\t',$_);
            my $line_A_size = @line_A;
           

            # test for the correct headers and extract chr id and chr length
            # load three columns of data into huge arrays
	    
				
	    if ($line_A[0] eq '@HD'){
		print "Found the SAM header and the following chromosome IDs and lengths: \n";
					}
			
	    elsif ($line_A[0] eq '@SQ'){
		@SAM_SN = split(':',$line_A[1]);
		@SAM_LN = split(':',$line_A[2]);
		push (@SAM_chr_id, $SAM_SN[1]);
		push (@SAM_chr_length, $SAM_LN[1]);
		print "Chromosome ID: $SAM_SN[1], Length: $SAM_LN[1] bp \n";
			}
	    elsif ($line_A[0] eq '@PG'){
		print "End of the SAM header.\n";
					}
	    #test for a forward strand alignment and extract 5' position
	    elsif ($line_A_size >3 && $line_A[1] == 0){
		push (@test_ID,$line_A[2]);
		push (@test_pos,$line_A[3]);
		
		$read_counter ++;
				}
	    #test for a reverse strand alignment and extract 3' position by adding the read length (STUPID SAM format!!!)		
	    elsif ($line_A_size >3 && $line_A[1] == 16){
		push (@test_ID,$line_A[2]);
		push (@test_pos,$line_A[3]+$read_length);
		
		$read_counter ++;
				}

        }

	# close in file handle
        close(IN);
	
	my $chr_list_size = @SAM_chr_id;
	my $data_list_size = @test_ID;
	print "\n Extracted data from $read_counter reads. \n";
	
		
#######################################################################################
# BUILD AN ARRAY MAP
#######################################################################################

print "\n\nIndexing all the data according to chromosome ID\n";
my $map_count = 0; # a counter variable

# Set bottom
$SAM_map{$test_ID[$map_count]} = 0;

$map_count ++;

# scan through the @test_ID array and mark the places where each new chromsomome starts

until ($map_count == $data_list_size){
  
      if ($test_ID[$map_count] ne $test_ID[$map_count-1]){
      
      $SAM_map{$test_ID[$map_count]} = $map_count;
      $map_count ++;
      
      }
      else{
      
	  $map_count ++;
	  
	  }

}
# output the number of chromosome types found as the number of hash keys.
$SAM_data_count = keys %SAM_map;
$mapsize = $map_count;

print "The data contained values corresponding to: $SAM_data_count chromosome(s)\n\n";



################################################################################
# Plot histogram of all read start positions for each chromosome
################################################################################

my $chr_counter = 0; # a counter variable
$map_count =0; #reset $map_count


until ($chr_counter == $chr_list_size){
  
	# set top bin of histogram
	my $top = $SAM_chr_length[$chr_counter];
	
	# define outfile name from infile name
        my $outfile = $SAM_chr_id[$chr_counter]."_".substr($infile_A,0,-4)."_dIEL_".$bin_width;
        $outfile .= '.sgr';
	
        
        # define new array to store required  read end position values
        my @read_pos;
        my $probe = 0;

        
        # use bin map to retrieve values
        
        $map_count = $SAM_map{$SAM_chr_id[$chr_counter]}; 
        
        while($SAM_chr_id[$chr_counter] eq $test_ID[$map_count]){
	    #test for end of data
            if ($map_count>=$mapsize -1){
            last;
            }      
                       
		push(@read_pos,$test_pos[$map_count]);
		$probe ++;
		$map_count ++;
		   

				

        
   }
 

        # Tally counter to plot histogram
		
		my $readarray_size= @read_pos;
		
		# Define the number of bins for the relevant chromosome
		my $bin_no = (int($top/$bin_width))+1;
		
		# Define the distribution frequency array
		my @dist_freq;
		my $i=0;
		
		# Fill the frequency distribution "bins" with zeroes
		for ($i=0; $i<$bin_no; $i++){
			push (@dist_freq, 0);
			}
			
		# Reset the incrementor and define the bin hit variable
		$i=0;
		my $bin_hit = 0;
		
		# The tally counter 
		while ($i < $readarray_size){
			$bin_hit = int($read_pos[$i]/$bin_width);
			$dist_freq[$bin_hit] ++;
			$i ++;
			}
		
		# Calculate the 3 bin moving average
		my @moving_ave;
		my $ma = 0;
		my $count = 1;
		push (@moving_ave,0);
		
		while ($count<$bin_no-1){
			$ma = (($dist_freq[$count--] + $dist_freq[$count] + $dist_freq[$count++])/3);
			push (@moving_ave,$ma);
			$count ++;
			}
			push (@moving_ave,0);
			
		
		
        # try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
        
            # print required data in tab-delimited format to output file
            # NK modifies to output chrn, bin and ma only
			for ($i=0; $i<$bin_no; $i++){
			
            print(OUT $SAM_chr_id[$chr_counter]."\t".($i*$bin_width)."\t".round($moving_ave[$i])."\n");
			}
            $chr_counter ++;
            print "Output: $outfile having found $probe data points.\n";
        }
        
        
         

   # close out file handle
        close(OUT);
        
       
}



}


print "end:\t", `date`."\n";
