#!/usr/bin/perl
# Written: Nick Kent, 12th Nov 2012
# Last updated: Nick Kent, 25th April 2013
# Last updated: Nick Kent, 23rd Oct 2019 to use Cwd
# USAGE:- perl Spombe_Chr_changer.plx
#
# This script takes chr1, chr2 format .sgr files that display in IGB using the
# S_pombe_Sep_2007 Genome Version, and converts the chromosome format to I, II etc allowing
# you to display the .sgrs in IGB using the S_pombe_May_2012 Genome Version. 
# 
################################################################################

use strict;
use warnings;
use Cwd;


################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $sgr_indir_path   - The directory containing the .sgr file
# $outdir_path  - The directory to store the flipped .sgr output file
# $sgr_outfile_name - User-defined name for the .sgr output file
# %chr_pairs - a hash containing the chr conversion criteria
#
################################################################################

my $sgr_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my %chr_pairs = (
'chr1'=>'I',
'chr2'=>'II',
'chr3'=>'III',
'chr4'=>'MT',
'chr5'=>'MTR');

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables and hash for the numeral swap


my $infile_sgr;
my $sgr_outfile;
my @line_sgr;
my @files_sgr;
my $sgr_size;


################################################################################
# Read in the .sgr file values to three enormous arrays
################################################################################


opendir(DIR,$sgr_indir_path) || die "Unable to access file at: $sgr_indir_path $!\n";

@files_sgr = readdir(DIR);

# process the input file within sgr_indir_path
foreach $infile_sgr (@files_sgr){    

    # ignore hidden files and only get those ending .sgr
    if (($infile_sgr !~ /^\.+/) && ($infile_sgr =~ /.*\.sgr/)){
        
       # define outfile name from infile name
        $sgr_outfile = substr($infile_sgr,0,-4)."_roman";
        $sgr_outfile .= '.sgr';

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

            # store the columns we want in the three new arrays and rename the chr column
			

            push(@sgr_chr,$chr_pairs{$line_sgr[0]});         
            push(@sgr_bin,$line_sgr[1]);
			push(@sgr_freq,($line_sgr[2]));
        }
        
        # close in file handle
        close(IN);
	

# store size of bin array
        $sgr_size = @sgr_freq;

print "Contains a whopping: $sgr_size bin values\n";


######################################################################################
# The output files
######################################################################################


# try and open the .sgr output file
        open(OUT,"> $outdir_path/$sgr_outfile")
             || die "Unable to open $sgr_outfile: $!";
        
print "Have just created $sgr_outfile\n";

# a counter variables
my $count = 0; # Counter for each line ID

until ($count == $sgr_size){ #until 1

print(OUT 
		$sgr_chr[$count]."\t".
		$sgr_bin[$count]."\t".
		$sgr_freq[$count]."\n");
		
		$count++;

} #until 1 closer


 # close .sgr out file handle
        close(OUT);


}}
