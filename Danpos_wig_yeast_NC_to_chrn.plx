#!/usr/bin/perl
# Written: Nick Kent, 10th Apr 2021
#
# USAGE:- perl Danpos_wig_yeast_NC_to_chrn.plx
#
# This script takes NC format .wig files from Danpos 2.2.2 processing (dpos) and
# changes chromosome ids to chr1, chr2 etc. format for
# Hughes 2008 build of the yeast genome.
#
# WARNING - the .wig format here will parse using IGB 9+ but not IGB 6; Not tested
# with UCSC or IGV - you may need to add a Header text to the output file.
# 
################################################################################

use strict;
use warnings;
use Cwd;


################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $wig_indir_path   - The directory containing the .wig file
# $outdir_path  - The directory to store the flipped .wig output file
# %chr_pairs - a hash containing the chr conversion criteria - check it's right!!!
################################################################################

my $wig_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my %chr_pairs = (
'chrom=gi|330443391|ref|NC_001133.9|'=>'chrom=chr1',
'chrom=gi|330443482|ref|NC_001134.8|'=>'chrom=chr2',
'chrom=gi|330443489|ref|NC_001135.5|'=>'chrom=chr3',
'chrom=gi|330443520|ref|NC_001136.10|'=>'chrom=chr4',
'chrom=gi|330443531|ref|NC_001137.3|'=>'chrom=chr5',
'chrom=gi|330443543|ref|NC_001138.5|'=>'chrom=chr6',
'chrom=gi|330443578|ref|NC_001139.9|'=>'chrom=chr7',
'chrom=gi|330443590|ref|NC_001140.6|'=>'chrom=chr8',
'chrom=gi|330443595|ref|NC_001141.2|'=>'chrom=chr9',
'chrom=gi|330443638|ref|NC_001142.9|'=>'chrom=chr10',
'chrom=gi|330443667|ref|NC_001143.9|'=>'chrom=chr11',
'chrom=gi|330443681|ref|NC_001144.5|'=>'chrom=chr12',
'chrom=gi|330443688|ref|NC_001145.3|'=>'chrom=chr13',
'chrom=gi|330443715|ref|NC_001146.8|'=>'chrom=chr14',
'chrom=gi|330443743|ref|NC_001147.6|'=>'chrom=chr15',
'chrom=gi|330443753|ref|NC_001148.4|'=>'chrom=chr16',
'chrom=gi|6226515|ref|NC_001224.1|'=>'chrom=chr17',
'chrom=gi|11466067|ref|NC_001398.1|'=>'chrom=chr18');

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables and hash for the numeral swap


my $infile_wig;
my $wig_outfile;
my @line_wig;
my @files_wig;
my $wig_size;


################################################################################
# Read in the .wig file values to three enormous arrays
################################################################################


opendir(DIR,$wig_indir_path) || die "Unable to access file at: $wig_indir_path $!\n";

@files_wig = readdir(DIR);

# process the input file within wig_indir_path
foreach $infile_wig (@files_wig){    

    # ignore hidden files and only get those ending .wig
    if (($infile_wig !~ /^\.+/) && ($infile_wig =~ /.*\.wig/)){
        
       # define outfile name from infile name
        $wig_outfile = substr($infile_wig,0,-4)."_chrn";
        $wig_outfile .= '.wig';

print "Found, and processing, $infile_wig \n";

open(IN, "$wig_indir_path/$infile_wig")
            || die "Unable to open $infile_wig: $!";
        
        # define three new arrays to store the .wig values from infile
        my @wig_out;
        my $new_fS;
        
        # loop through infile to get values
        while(<IN>){

            chomp;
            
            if($_ !~ m/fixed/){
            push (@wig_out,$_);
            }
            else{
            
            # split line by delimiter and store elements in an array
            @line_wig = split('\s',$_);
            
            # brute find and replace using hash
            
            $new_fS ="$line_wig[0] $chr_pairs{$line_wig[1]} $line_wig[2] $line_wig[4] $line_wig[5]";
            push (@wig_out,$new_fS);
            }
           
        }
        
        # close in file handle
        close(IN);
	

# store size of bin array
        $wig_size = @wig_out;

print "Contains a whopping: $wig_size bin values\n";


######################################################################################
# The output files
######################################################################################


# try and open the .wig output file
        open(OUT,"> $outdir_path/$wig_outfile")
             || die "Unable to open $wig_outfile: $!";
        
print "Have just created $wig_outfile\n";

# a counter variables
my $count = 0; # Counter for each line ID

until ($count == $wig_size){ #until 1

print(OUT 
		$wig_out[$count]."\n");
		
		$count++;

} #until 1 closer


 # close .wig out file handle
        close(OUT);


}}
