#!/usr/bin/perl
# Written: Nick Kent, 12th Sept 2010
# Modified from Site_Writer_CFD. plx: Nick Kent 19th Apr 2012 - solve over-run error.
# Modified: Nick Kent 20th Nov 2018 - force numeric to string conversion to allow numeric Chr IDs
# Modified: Nick Kent 22nd Feb 2019 - output file to 2 col raw text; input file ends changed; .bedGraph input
# Modified: Nick Kent 22nd May 2019 - revert to .sgr and .txt input; change to ignore non-matched "chr" site records
# from site wrangler processing.

# USAGE:- perl CFD_Plotter.plx
#
# NOTE Nick needs to error trap the non-matching chr ID thing
#
# FUNCTION:
# This script takes .txt files containing a list of sites/genomic features (these
# could be TSSs or TF sites or whatever you want) and compares it with 
# whole-genome, Partn .sgr files. It then outputs CUMULATIVE FREQUENCY DISTRIBUTION
# values over a user-specified bin range centered on, and surrounding the sites.
#
# The output file can be used to plot average chromatin particle environments for 
# different sorts of TSS for example. Sites close to chromosome ends, which would 
# not yield the full range of data are ignored, but reported at the command line.
#
#
# INPUT AND OUTPUT (all tab-delimited):
# The input .txt files should have four columns: chrn; Site ID; site pos; strand.
#
# The input .sgr files should have three columns: chrn; bin pos; paired-read dyad freq.
#
# The output .txt file has an input file header and column headers and returns 
# 2 columns: Bin (relative to Site); normalised cumulative freq. The idea is to plot the 
# columns as a line graph to produce a TREND GRAPH for the data. Each bins F+R
# frequencies are normalised to the average F+R frequency for the entire bin window.
# 
# Note: Use multiple .sgrs and then cat or paste the CFD.txt files for processing in R or 
# Excel - particularly useful for plotting surface landscape graphs.
#
# Note: The script handles F and R strand data separately. If you give it all F (or
# all R) strand sites it will work just fine, however, it will also throw a load of
# uninitialised variable warnings at the command line. If you find this upsetting,
# stick a # in front of use warnings (below)
# 
# For development see: Kent et al.,(2011) Chromatin particle spectrum analysis: a 
# method for comparative chromatin structure analysis using paired-end mode 
# next-generation DNA sequencing. NAR 39: e26.
################################################################################

use strict;
use warnings;
use Cwd;
use List::Util;

################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $sgr_indir_path   - The directory containing the full genome Partn.sgr files
# $siteID_indir_path - The directory containing the site list .txt file
# $outdir_path  - The directory to store the output files
#
# $bin_window - number of bins surrounding the site of interest. E.g. if you set
# this to 40 then you will get 40 bins either side of your site - 400bp if you
# were using 10bp binned data.
#
# $bin_size - binning interval of .sgr file in base pairs.
#
# $output_scale - controls how many bins are included in the output file. If set
# to 1 you will get every bin (use this). Set to 3 to output only every third bin
# in the series.You can use this feature to scale output files derived from input 
# .sgr data with different bin intervals.
################################################################################

my $sgr_indir_path =cwd."/Sgr_in";
my $siteID_indir_path =cwd."/Site_in";
my $outdir_path =cwd."/out";
my $bin_window = 120;
my $bin_size = 10;
my $output_scale = 1;

################################################################################
################################################################################
# MAIN PROGRAM 
################################################################################
################################################################################

# define some variables
my $infile_sgr;
my $infile_siteID;
my $cfd_outfile;
my @line_siteID; 
my @line_sgr;
my @files_sgr;
my @files_siteID;
my $sgr_size;
my $F_siteID_size;
my $R_siteID_size;
my %bin_map; 
my $chr_count;
my $descriptor;
my $spread = $bin_window*$bin_size;

print "\nScript is set to calculate normalised cumulative sequence read mid-point frequency values: $spread bp either side of the genomic locations/sites you have provided; binned every $bin_size bp.\n\n";
################################################################################
# Get site list and write to an array - from .txt format with four columns:
# chrn;siteID; site position; F/R 
################################################################################

# store input file name in an array
opendir(DIR,$siteID_indir_path) || die "Unable to access file at: $siteID_indir_path $!\n";

@files_siteID = readdir(DIR);

# process the input file within siteID_indir_path
foreach $infile_siteID (@files_siteID){    

    # ignore hidden files and only get those ending .txt
    if (($infile_siteID !~ /^\.+/) && ($infile_siteID =~ /.*\.txt/)){
        
    $descriptor = substr($infile_siteID,0, -4);

        
print "Found, and processing, $infile_siteID \n";

open(IN, "$siteID_indir_path/$infile_siteID")
            || die "Unable to open $infile_siteID: $!";
        
        # define strand-specific arrays to store site chromosome no., and position 
		my @F_site_chr;
		my @F_site_pos;
		my @R_site_chr;
		my @R_site_pos;
	
	# loop through infile to get values
        while(<IN>){

	    chomp;

	    # split line by delimiter and store elements in an array
            @line_siteID = split('\t',$_);

            # store the required chrn, position in two pairs of strand-specific
			# arrays
			# test for unmatch "chr" data and skip
			
			unless($line_siteID[0] eq "chr"){
			
			if($line_siteID[3] =~ "F"){ # infile if 1
			#altered here to force numeric/string
				push(@F_site_chr,"".$line_siteID[0]);
				push(@F_site_pos,$line_siteID[2]);
				
			}elsif($line_siteID[3] =~ "R"){
			#altered here to force numeric/string
				push(@R_site_chr,"".$line_siteID[0]);
				push(@R_site_pos,$line_siteID[2]);
				
			}else{
			
				print "Failed to match strand at $line_siteID[0], $line_siteID[1], $line_siteID[2]\n";
			
			  } #infile if 1 closer
            } #unless closer
        }

	# close in file handle
        close(IN);
	closedir(DIR);
        
        # store sizes of the arrays
        $F_siteID_size = @F_site_pos;
		$R_siteID_size = @R_site_pos;


print "Contains: $F_siteID_size forward strand site IDs; $R_siteID_size reverse
 strand site IDs\n";


################################################################################
# Read in the .sgr file values to three enormous arrays
################################################################################


opendir(DIR,$sgr_indir_path) || die "Unable to access file at: $sgr_indir_path $!\n";

@files_sgr = readdir(DIR);

# process the input file within sgr_indir_path
foreach $infile_sgr (@files_sgr){   

# define some arrays that will be reset during each iteration
my @F_cfd_bin;
my @F_cfd_freqsum;
my @R_cfd_bin;
my @R_cfd_freqsum; 

    # ignore hidden files and only get those ending .sgr
    if (($infile_sgr !~ /^\.+/) && ($infile_sgr =~ /.*\.sgr/)){
        
       

print "Found, and processing, $infile_sgr \n";

open(IN, "$sgr_indir_path/$infile_sgr")
            || die "Unable to open $infile_sgr: $!";
        
        # define three new arrays to store the .sgr values from infile
        my @sgr_chr;
        my @sgr_bin;
		my @sgr_freq;
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_sgr = split('\t',$_);

            # store the columns we want in the three new arrays
            #altered here to force numeric/string
            push(@sgr_chr,"".$line_sgr[0]);
            push(@sgr_bin,$line_sgr[1]);
			push(@sgr_freq,$line_sgr[2]);
        }
        
        # close in file handle
        close(IN);
	

# store size of bin array
        $sgr_size = @sgr_freq;

print "Contains a whopping: $sgr_size bin values\n";


#######################################################################################
# BUILD THE BIN MAP
#######################################################################################

my $map_count = 0; # a counter variable

# Set bottom
$bin_map{$sgr_chr[$map_count]} = 0;

$map_count ++;

# scan through the @sgr_chr array and mark the bins where each new chromsomome starts

until ($map_count == $sgr_size){
  
      if ($sgr_chr[$map_count] ne $sgr_chr[$map_count-1]){
      
      $bin_map{$sgr_chr[$map_count]} = $map_count;
      $map_count ++;
      
      }
      else{
      
	  $map_count ++;
	  
	  }

}
# output the number of chromosome types found as the number of hash keys.
$chr_count = keys %bin_map;

print "Contains values for: $chr_count chromosomes\n";

#######################################################################################
# FORWARD STRAND .sgr calculations:

# some counter variables
my $site_count = 0; # Counter for each site ID
my $bin_count = 0; # Counter .sgr bin numbers
my $cfd_count = 0; # Counter for the cfd arrays
my $top_limit = 0; # A top limit for $bin_window
my @F_out_chr; # F .sgr output array for chr
my @F_out_bin; # F .sgr output array for bin pos
my @F_out_freq; # F .sgr output array for read freq
my $F_out_size = 0; # Size of F .sgr output arrays
my $i=0; # An iterator variable

until ($site_count == $F_siteID_size){ #until 1

# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($F_site_pos[$site_count]/$bin_size) + $bin_map{$F_site_chr[$site_count]}) - 3;

# this looks mad, but it allows me to recycle all the code from the last version, and takes up any
# rounding slack which would come from different $bin_size values




# find an .sgr bin which contains the current site
	until ($F_site_chr[$site_count] eq $sgr_chr[$bin_count] &&
	       $F_site_pos[$site_count] >= $sgr_bin[$bin_count] &&
	       $F_site_pos[$site_count] <  $sgr_bin[$bin_count +1]){ #until 2

	       
			$bin_count ++;
			


	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($F_site_chr[$site_count] ne $sgr_chr[$bin_count] ||
	   $F_site_chr[$site_count] ne $sgr_chr[$top_limit]){ #if 1

	print "Can't output forward strand values for $F_site_chr[$site_count] site: $F_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the F .sgr arrays and add values to  F cfd freq array

	until ($bin_count == $top_limit+1){ #until 3

		push (@F_out_chr,$sgr_chr[$bin_count]);
		push (@F_out_bin,$sgr_bin[$bin_count]);
		push (@F_out_freq,$sgr_freq[$bin_count]);

		$F_cfd_freqsum[$cfd_count] += $sgr_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$F_out_size = @F_out_chr;


		
####################################################################################		
# REVERSE STRAND .sgr calculations:

# reset the counter variables and define some more arrays
$site_count = 0; # Counter for each site ID
$cfd_count = 0; # Counter for the cfd arrays
$bin_count = 0;
my @R_out_chr; # R .sgr output array for chr
my @R_out_bin; # R .sgr output array for bin pos
my @R_out_freq; # R .sgr output array for read freq
my $R_out_size = 0; # Size of F .sgr output arrays

until ($site_count == $R_siteID_size){ #until 1

# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($R_site_pos[$site_count]/$bin_size) + $bin_map{$R_site_chr[$site_count]}) - 3;


# find an .sgr bin which contains the current site
	until ($R_site_chr[$site_count] eq $sgr_chr[$bin_count] &&
	       $R_site_pos[$site_count] >= $sgr_bin[$bin_count] &&
	       $R_site_pos[$site_count] <  $sgr_bin[$bin_count +1]){ #until 2


			
			$bin_count ++;
			


	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($R_site_chr[$site_count] ne $sgr_chr[$bin_count] ||
	   $R_site_chr[$site_count] ne $sgr_chr[$top_limit]){ #if 1

	print "Can't output reverse strand values for $R_site_chr[$site_count] site: $R_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the R .sgr arrays and add values to  R cfd freq array


	until ($bin_count == $top_limit+1){ #until 3

		push (@R_out_chr,$sgr_chr[$bin_count]);
		push (@R_out_bin,$sgr_bin[$bin_count]);
		push (@R_out_freq,$sgr_freq[$bin_count]);

		$R_cfd_freqsum[$cfd_count] += $sgr_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$R_out_size = @R_out_chr;



######################################################################################
# The output file 
######################################################################################

# define outfile name and set correct endings


	$cfd_outfile = substr($infile_sgr,0,-4)."_".$descriptor."_CFD";
       $cfd_outfile .= '.txt';

# try and open the .cfd output file
        open(OUT,"> $outdir_path/$cfd_outfile")
             || die "Unable to open $cfd_outfile: $!";
        


# Set counter variables and define new arrays
$bin_count = 0;
$cfd_count = 0;
my $cfd_sum = 0; # a sum of sums for normalizing the data
my $norm_factor = 0; # calced from $cfd_sum
my $R_cfd_count = $bin_window*2;
my @FandR_cfd; # array to hold summed F and R strand CFD values
my @R_cfd; # array to hold ordered R strand CFD values
my @normalized_cfd;

$bin_count -= $bin_window;

	until ($bin_count == $bin_window+1){ #until 4

		# re-order reverse strand cfd freqsum values
		push (@R_cfd, $R_cfd_freqsum[$R_cfd_count]); 
		
		# calculate summed value for both F and R cfd freqsums
		push (@FandR_cfd, $F_cfd_freqsum[$cfd_count] + $R_cfd_freqsum[$R_cfd_count]);
	
		$bin_count ++;
		$cfd_count ++;
		$R_cfd_count --;

		} #until 4 closer

		
# Need to find average read values over bin_window to normalize data

$cfd_sum += $_ for @FandR_cfd;
$norm_factor = $cfd_sum/(($bin_window*2)+1);

# reset counters once more
$bin_count = (0-$bin_window);
$cfd_count = 0;

# print a header for the CFD.txt file so you can read it in Excel

print (OUT "Values from $cfd_outfile\n");
print (OUT "CFD sum: $cfd_sum\n");
print (OUT "Normalization Factor: $norm_factor\n");

# print column headers

print (OUT "Bin (bp)"."\t"."Normalised Cumulative Read Freq"."\n");

# print data values
		until ($bin_count >= $bin_window+1){ #until 5

		print(OUT 
		$bin_count*$bin_size."\t".
		$FandR_cfd[$cfd_count]/$norm_factor."\n");
	
	
		$bin_count += $output_scale;
		$cfd_count += $output_scale;

		} #until 5 closer
		
 # close .cfd out file handle
        close(OUT);

        print "Have just created $cfd_outfile\n\n";

}}}
}



print "Finished.\n\n";
