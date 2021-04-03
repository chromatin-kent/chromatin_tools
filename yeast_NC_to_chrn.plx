#!/usr/bin/perl
# Written: Nick Kent, 12th Nov 2012
# Last updated: Nick Kent, 12th Nov 2012

# USAGE:- perl yeast_toSacCer3.plx
#
# This script takes NC format .sgr files from SAM2PartN_sgr and converts to chr1, chr2 format for
# Hughes 2008 build of the yeast genome.
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
# %chr_pairs - a hash containing the chr conversion criteria - check it's right!!!
#
################################################################################

my $sgr_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my %chr_pairs = (
'gi|330443391|ref|NC_001133.9|'=>'chr1',
'gi|330443482|ref|NC_001134.8|'=>'chr2',
'gi|330443489|ref|NC_001135.5|'=>'chr3',
'gi|330443520|ref|NC_001136.10|'=>'chr4',
'gi|330443531|ref|NC_001137.3|'=>'chr5',
'gi|330443543|ref|NC_001138.5|'=>'chr6',
'gi|330443578|ref|NC_001139.9|'=>'chr7',
'gi|330443590|ref|NC_001140.6|'=>'chr8',
'gi|330443595|ref|NC_001141.2|'=>'chr9',
'gi|330443638|ref|NC_001142.9|'=>'chr10',
'gi|330443667|ref|NC_001143.9|'=>'chr11',
'gi|330443681|ref|NC_001144.5|'=>'chr12',
'gi|330443688|ref|NC_001145.3|'=>'chr13',
'gi|330443715|ref|NC_001146.8|'=>'chr14',
'gi|330443743|ref|NC_001147.6|'=>'chr15',
'gi|330443753|ref|NC_001148.4|'=>'chr16',
'gi|6226515|ref|NC_001224.1|'=>'chr17',
'gi|11466067|ref|NC_001398.1|'=>'chr18');

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
        $sgr_outfile = substr($infile_sgr,0,-4)."_chrn";
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
