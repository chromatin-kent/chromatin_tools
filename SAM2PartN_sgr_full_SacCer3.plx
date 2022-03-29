#!/usr/bin/perl
use strict;
use warnings;
use Math::Round;
use Cwd;
#################################################################################
# Written by Nick Kent, Dec 2017
# Fiddled with: Nick Kent, Dec 20th 2017 to add more processing information and time stamps.
# Fiddled with: Nick Kent, Dec 21st 2017 to load SAM data into memory.
# Fiddled with: Nick Kent, Jan 4th 2018 to process multiple PartN values - seems to work
# Fiddled with: Nick Kent, Jan 5th 2018 to solve end-of-data run-off error
# Fiddled with: Nick Kent, Oct 23rd 2019 to use cwd
# Fiddled with: Nick Kent, Mar 29th 2022 to run under Gomphus environment for BI3008
#################################################################################
# USAGE:- perl SAM2PartN_sgr_full.plx
#
# This script will take a SORTED* Bowtie 1 CPSA paired read alignment .sam format file
# and will generate chromosome-specific .sgr files containing 3MA-smoothed read mid-
# point frequency values for MNase-resistant chromatin particle reads.
#
# This script will NOT run on a weedy Windows Laptop, but requires a high memory
# Linux machine. Typically you will need RAM eqivalent to the 0.75 * size 
# of your .sam file. There is also a limit to the number of reads that can be processed,
# and the number of genomic bins. These limits are both = 2^31 - the max size of a perl array. This
# script was designed for model eukaryote (yeast, Drosophila, Arabidopsis) and prokaryotic
# genomes, but may not work for mammalian size genomes.
#
# User specifies a list of  "PartN" value in bp - e.g. 150 for a nucleosome, 50 for a TF.
# User specifies a window +/- PartN as a fraction of 1 (i.e. 0.2 is eqivalent to +/- 20%)
#
# SUGGESTION: once finished the user could "cat" (bash shell) the individual .sgr 
# files into one "whole genome" file, and use namechanger.plx to alter Chromosome IDs.
#
# * To sort use: ./samtools sort Input.sam -o Input_sorted.sam
#
# NOTE: this script is the bastard love-child of pair_read_histo.plx, histogram.plx and SiteWriterCFD.plx
#
# To do: MORE ERROR TRAPPING (e.g. non-sorted .sam; general fuck-wittery;); rename test variables.
################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $indir_path   - The directory containing the SORTED .sam file to be processed
# $outdir_path  - The directory to store the .sgr output files
# $bin_width    - The histogram bin width value in bp (10bp is good for most things)
# @partn 	- A bracketed, comma delimited list of "chromatin particle sizes" e.g. (50,150,300)
# $pwind	- The "particle window" e.g. 0.2 = +/- 20% of $partn
################################################################################

my $inA_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my $bin_width = 10;
my @partn = (50,150,300);
my $pwind = 0.2;

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
my @test_ISIZE;
my $read_counter = 0;
my %SAM_map; 
my $SAM_data_count;
my $partn_size = @partn;
my $partsize;
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
print "\nThe following MNase protection particle sizes have been specified:\n ";
    
    foreach $partsize (@partn){
	print "$partsize bp\n ";
    
			}
print "\nParticle size window will be +/-".($pwind*100)." percent.\n";
print "\nFound, and processing, $infile_A \n\n";


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
        print "Starting to read SAM data into memory - this might take a while; please be patient.\n";
					}
			
	    elsif ($line_A_size >3 && $line_A[8]>0){
		push (@test_ID,$line_A[2]);
		push (@test_pos,$line_A[3]);
		push (@test_ISIZE,$line_A[8]);
		$read_counter ++;
				}
				
	  elsif ($line_A_size >3 && $line_A[8]<=0){
				}
				
	  else{
		print "\n Not sure this is a normal SAM file. Will stop so that you can re-check \n";
		exit;
	  }
        }

	# close in file handle
        close(IN);
	
	my $chr_list_size = @SAM_chr_id;
	my $data_list_size = @test_ID;
	
	
		
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
# Cycle through each PartN size
################################################################################
foreach $partsize (@partn){

# define outfile name from infile name NOTE this changes original multiple file output to concatenation
my $outfile = substr($infile_A,0,-4)."_Part".$partsize."_".$bin_width;
$outfile .= '.sgr';

# try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        


# define lowwer and upper values of ISIZE as a percentage from $pwind
my $ISIZE_low = ($partsize - ($pwind * $partsize));
my $ISIZE_high =($partsize + ($pwind * $partsize));


################################################################################
# Plot histogram of PartN midpoint position for each chromosome
################################################################################

my $chr_counter = 0; # a counter variable
$map_count =0; #reset $map_count


until ($chr_counter == $chr_list_size){
  
	# set top bin of histogram
	my $top = $SAM_chr_length[$chr_counter];
	
        
        # define new array to store required dyad (paired read mid-point) position values
        my @dyad_pos;
        my $probe = 0;

        
        # use bin map to retrieve values
        
        $map_count = $SAM_map{$SAM_chr_id[$chr_counter]}; 
        
        while($SAM_chr_id[$chr_counter] eq $test_ID[$map_count]){
	    #test for end of data
            if ($map_count>=$mapsize -1){
            last;
            }      
	    #test for match within pwindow
            if ($test_ISIZE[$map_count] > $ISIZE_low && $test_ISIZE[$map_count] < $ISIZE_high){
                       
		push(@dyad_pos,($test_pos[$map_count] + ($test_ISIZE[$map_count] * 0.5)));
		$probe ++;
		$map_count ++;
		   
	      }
	  else{
		$map_count ++;
		}
				

        
   }
 

        # Tally counter to plot histogram
		
		my $dyadarray_size= @dyad_pos;
		
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
		while ($i < $dyadarray_size){
			$bin_hit = int($dyad_pos[$i]/$bin_width);
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
        
            # print required data in tab-delimited format to output file
            # NK modifies to output chrn, bin and ma only
			for ($i=0; $i<$bin_no; $i++){
			
            print(OUT $SAM_chr_id[$chr_counter]."\t".($i*$bin_width)."\t".round($moving_ave[$i])."\n");
			}
            $chr_counter ++;
            print "Output: $outfile having found $probe data points for $partsize bp particles.\n";
        }
        
        
         
 }

        
       
}

   # close out file handle
        close(OUT);

}


print "end:\t", `date`."\n";
