#!/usr/bin/perl
# Written: Nick Kent, 16th June 2011
# Last updated 16th June 2011
# Last updated 5th Sep 2019 cwd
#
# USAGE:- perl interval_wrangler_ID.pl
#
# This script is a modifcation of interval_comparison.pl. It takes an interval file
#  i.e. introns or intergenic regions whatever and asks whether bins in an .sgr file
# output by Peak_Comparison.pl occur within that region. 
#
# The script takes two input files in tab-delimited .txt format:
#
# The peak_file is a standard .sgr: chrn; summit(peak) bin; summit bin read frequency.
#
# The interval file should contain columns: intervalID no.; 
# interval descriptor (ID); chrn; start pos; end pos.
# e.g.: intron3; YALintron; chr1; 2346; 2789.
#
# The Script outputs TWO files:
# (i) a three column .sgr file for values present within the intervals
# (ii) a four column .txt file which shows ID, interval ID, chrn, summit bin
# It also outputs at command line the total number of peaks which found a matching interval. 
#
################################################################################

use strict;
use warnings;
use Cwd;

################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $interval_indir_path   - The directory containing the interval .txt file
# $peak_indir_path   - The directory containing the peakmarker .txt file
# $outdir_path  - The directory to store the .txt output file
# $no - The values entered to output file if no peak is matched to an interval
# $interval_outfile_name - An informative output file name for the output .txt
################################################################################

my $interval_indir_path =cwd."/interval_file";
my $peak_indir_path =cwd."/peak_file";
my $outdir_path =cwd."/out";
my $cwd = getcwd;
my $no = 0;
my $interval_outfile_name = "WT_NOT_chd1_MalabatTSS-350_interval_Part50";
    
################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables

my $infile_interval;
my $infile_peak;
my $outfile; 
my @line_interval;
my @line_peak;
my @files_interval;
my @files_peak;
my $interval_size;
my $peak_size;

################################################################################
# Read in the interval file values to five arrays
################################################################################

# store input file name in an array
opendir(DIR,$interval_indir_path) || die "Unable to access file at: $interval_indir_path $!\n";

@files_interval = readdir(DIR);

# process the input file within indir_path
foreach $infile_interval (@files_interval){    

    # ignore hidden files and only get those ending .txt
    if (($infile_interval !~ /^\.+/) && ($infile_interval =~ /.*\.txt/)){
        
        
print "Found, and processing, $infile_interval \n";

open(IN, "$interval_indir_path/$infile_interval")
            || die "Unable to open $infile_interval: $!";
        
        # define the arrays to store required values from infile
        my @interval_chr; # the interval chrn
	my @interval_id; # an informative interval ID number
        my @interval_start; # the start position of the interval
	my @interval_end; # the end position of the interval
	my @interval_name; # A text descriptor if available/extractable
	
	
	
	# loop through infile to get values
        while(<IN>){
           
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_interval = split('\t',$_);

            # store the columns we want in five new arrays
            
	    push(@interval_id,$line_interval[0]);
	    push(@interval_name,$line_interval[1]);
	    push(@interval_chr,$line_interval[2]);
	    push(@interval_start,$line_interval[3]);
	    push(@interval_end,$line_interval[4]);
	    
            
        }

	# close in file handle
        close(IN);
	closedir(DIR);
        
        # store size of bin array
        $interval_size = @interval_chr;


print "Contains: $interval_size intervals\n";


################################################################################
# Read in the peak file values to four arrays
################################################################################


opendir(DIR,$peak_indir_path) || die "Unable to access file at: $peak_indir_path $!\n";

@files_peak = readdir(DIR);

# process the input file within indir_path
foreach $infile_peak (@files_peak){    

    # ignore hidden files and only get those ending .sgr
    if (($infile_peak !~ /^\.+/) && ($infile_peak =~ /.*\.sgr/)){
        
       

print "Found, and processing, $infile_peak \n";

open(IN, "$peak_indir_path/$infile_peak")
            || die "Unable to open $infile_peak: $!";
        
        # define three new arrays to store required values from .sgr infile
       
        my @peak_chr; # the .sgr file chrn
        my @peak_bins; # the peak summit bin in bp
        my @peak_freq; # the averaged read frequency reported in the peak summit bin
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_peak = split('\t',$_);

            # store the columns we want in the four new arrays
           
	    push(@peak_chr,$line_peak[0]);
            push(@peak_bins,$line_peak[1]);
            push(@peak_freq,$line_peak[2]);
        }
        
        # close in file handle
        close(IN);
	closedir(DIR);

# store size of bin array
        $peak_size = @peak_chr;

print "Contains: $peak_size peak values\n";


######################################################################################
# The .sgr output file
######################################################################################

# define outfile name and set ending to .txt
        $outfile = $interval_outfile_name;
        $outfile .= '.sgr';

# try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
print "Have just created $outfile\n";
        
# some counter variables
my $interval_count = 0; # interval site number
my $peak_count = 0; # Peak site number
my $hit_count = 0; # A counter for positive site/peak hits
my $match = 0; # running score of total peaks within this interval


# a very inefficient way of searching for matches - but it works
until ($interval_count == $interval_size){

	until ($peak_count == $peak_size){ 


	# this increments the peak counter if the chrn values do not match
	if ($interval_chr[$interval_count] ne $peak_chr[$peak_count]){

		$peak_count++;

        # this tests whether or not a peak bin occurs within an interval 

	}elsif ($peak_bins[$peak_count]>=($interval_start[$interval_count]) && 
			$peak_bins[$peak_count]<=($interval_end[$interval_count]) && 
			$interval_chr[$interval_count] eq $peak_chr[$peak_count]){

		# this writes out the 3 column tab-delimited .sgr
		print(OUT 
		
		$peak_chr[$peak_count]."\t".
		$peak_bins[$peak_count]."\t".
		$peak_freq[$peak_count]."\n"); # removed extra tab here
		
		$hit_count++;
		$peak_count++;

	}else{
		$peak_count ++;
	}
}
	
	if ($hit_count==0){

		$interval_count++;
		$peak_count = 0;

	}else{
		$match += $hit_count;
		$interval_count++;
		$hit_count = 0;
		$peak_count = 0;
	}


}
	
       
        # close out file handle
        close(OUT);


######################################################################################
# The .txt output file
######################################################################################

# define outfile name and set ending to .txt
        $outfile = $interval_outfile_name;
        $outfile .= '.txt';

# try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
print "Have just created $outfile\n";
        
# some counter variables
$interval_count = 0; # interval site number
$peak_count = 0; # Peak site number
$hit_count = 0; # A counter for positive site/peak hits
$match = 0; # running score of total peaks within this interval


# a very inefficient way of searching for matches - but it works...again
until ($interval_count == $interval_size){

	until ($peak_count == $peak_size){ 


	# this increments the peak counter if the chrn values do not match
	if ($interval_chr[$interval_count] ne $peak_chr[$peak_count]){

		$peak_count++;

        # this tests whether or not a peak bin occurs within an interval 

	}elsif ($peak_bins[$peak_count]>=($interval_start[$interval_count]) && 
			$peak_bins[$peak_count]<=($interval_end[$interval_count]) && 
			$interval_chr[$interval_count] eq $peak_chr[$peak_count]){

		# this writes out the 3 column tab-delimited .sgr
		print(OUT
 
		$interval_name[$interval_count]."\t".	
		$interval_id[$interval_count]."\t".	
		$peak_chr[$peak_count]."\t".
		$peak_bins[$peak_count]."\n");

		
		$hit_count++;
		$peak_count++;

	}else{
		$peak_count ++;
	}
}
	
	if ($hit_count==0){

		$interval_count++;
		$peak_count = 0;

	}else{
		$match += $hit_count;
		$interval_count++;
		$hit_count = 0;
		$peak_count = 0;
	}


}
	
       
        # close out file handle
        close(OUT);

	print "Out of $peak_size peaks, $match occured within the specified intervals\n";

}}}}
       
