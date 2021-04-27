#!/usr/bin/perl
# Written: Nick Kent, 12th Nov 2012
# Last updated: Nick Kent, 12th Nov 2012
# Modified: Nick Kent, 22nd Feb 2019
# Modified: Nick Kent, 27 Apr 2021
# USAGE:- perl Pombe_SITE_to_Roman.plx
#
# This script takes chr1, chr2 format site.txt files and converts to chrI, chrII format for
# the SacCer3 build of the yeast genome.
# 
################################################################################

use strict;
use warnings;
use Cwd;


################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $txt_indir_path   - The directory containing the .txt file
# $outdir_path  - The directory to store the flipped .txt output file
# $txt_outfile_name - User-defined name for the .txt output file
# %chr_pairs - a hash containing the chr conversion criteria - check it's right!!!
#
################################################################################

my $txt_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my %chr_pairs = (
'chr1'=>'I',
'chr2'=>'II',
'chr3'=>'III');

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables 


my $infile_txt;
my $txt_outfile;
my @line_txt;
my @files_txt;
my $txt_size;


################################################################################
# Read in the .txt file values to four enormous arrays
################################################################################


opendir(DIR,$txt_indir_path) || die "Unable to access file at: $txt_indir_path $!\n";

@files_txt = readdir(DIR);

# process the input file within txt_indir_path
foreach $infile_txt (@files_txt){    

    # ignore hidden files and only get those ending .txt
    if (($infile_txt !~ /^\.+/) && ($infile_txt =~ /.*\.txt/)){
        
       # define outfile name from infile name
        $txt_outfile = substr($infile_txt,0,-4)."_roman";
        $txt_outfile .= '.txt';

print "Found, and processing, $infile_txt \n";

open(IN, "$txt_indir_path/$infile_txt")
            || die "Unable to open $infile_txt: $!";
        
        # define  new arrays to store the .txt values from infile
        my @txt_chr;
        my @txt_siteid;
        my @txt_sitepos;
        my @txt_sitestrand;
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_txt = split('\t',$_);

            # store the columns we want in the  new arrays and rename the chr column
			

            push(@txt_chr,$chr_pairs{$line_txt[0]});         
            push(@txt_siteid,$line_txt[1]);
			push(@txt_sitepos,($line_txt[2]));
			push(@txt_sitestrand,($line_txt[3]));
        }
        
        # close in file handle
        close(IN);
	

# store size of bin array
        $txt_size = @txt_sitepos;

print "Contains a whopping: $txt_size bin values\n";


######################################################################################
# The output files
######################################################################################


# try and open the .txt output file
        open(OUT,"> $outdir_path/$txt_outfile")
             || die "Unable to open $txt_outfile: $!";
        
print "Have just created $txt_outfile\n";

# a counter variables
my $count = 0; # Counter for each line ID

until ($count == $txt_size){ #until 1

print(OUT 
		$txt_chr[$count]."\t".
		$txt_siteid[$count]."\t".
		$txt_sitepos[$count]."\t".
		$txt_sitestrand[$count]."\n");
		
		$count++;

} #until 1 closer


 # close .txt out file handle
        close(OUT);


}}
