#!/usr/bin/perl
# Written: Nick Kent, Nov 2011
# Last updated: Nick Kent, Nov 2011
# Last fiddled: Nick Kent, 23 Jan 2020 - use cwd
#
# Cluster 3.0 takes a C3.txt file from ClusterWriter_Full.plx and outputs
# k-means (or whatever) clusters of that data in ordered .cdt format for
# input into TreeView. Sometimes it would be nice if you could see the 
# exact same ordered dataset associated with a mutant condition. This would let
# you visualise changes in cluster composition between WT and mutant without
# reclustering and sorting the .kgg data. So we need a script which will take the 
# .cdt file from one experimental condition and create an identical .cdt file (in
# terms of ID_string ordering) BUT containing the data points from the C3.txt file
# corresponding to another experimental condition. This, is that script...
#
# Before you begin, make sure you have run CFD_Plotter_Full_C3_nonorm.plx using the SAME
# site_in file on both datasets - e.g. have done a WT vs mutant comparison to create 
# WT_C3.txt and mut_c3.txt output files. You would then use the WT_C3.txt file to 
# generate a WT.cdt file using Cluster 3.0. The WT.cdt file would then be used
# as a template to order the data from the mut_C3.txt file. Got that?
#
# The script takes two input files in tab-delimited .txt format:
#
# 1. The A file is the template (e.g. WT) .cdt file
#
# 2. The B file is the corresponding data to be reordered (e.g. mut) C3.txt file
#
#
# The script will output a .cdt format file with identical row and column headers
# to the A file but containing the B file data values
#
#  WARNING this is a development script (aren't they all...)
################################################################################

use strict;
use warnings;
use Cwd;


################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $inA_indir_path   - The directory containing the A .cdt file
# $inB_indir_path   - The directory containing the B _C3.txt file
# $outdir_path  - The directory to store the .cdt output file
#
################################################################################

my $inA_indir_path =cwd."/A_in";
my $inB_indir_path =cwd."/B_in";
my $outdir_path =cwd."/out";

    
################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables

my $infile_A;
my $infile_B;
my $outfile; 
my @line_A;
my @line_B;
my @files_A;
my @files_B;
my $size_A;
my $size_B;

################################################################################
# Read in the A file values to four arrays
################################################################################

# store input file name in an array
opendir(DIR,$inA_indir_path) || die "Unable to access file at: $inA_indir_path $!\n";

@files_A = readdir(DIR);

# process the input file within indir_path
foreach $infile_A (@files_A){    

    # ignore hidden files and only get those ending .cdt
    if (($infile_A !~ /^\.+/) && ($infile_A =~ /.*\.cdt/)){
        
        
print "Found, and processing, $infile_A \n";

open(IN, "$inA_indir_path/$infile_A")
            || die "Unable to open $infile_A: $!";
        
        # define the arrays to store required values from infile
       
		my @cdt_id; # the .cdt ORF ID strings
		my @cdt_lines; # the full .cdt lines
        
	
	# loop through infile to get values
        while(<IN>){
           
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_A = split('\t',$_);

            # store the columns we want in the two arrays

            push(@cdt_id,$line_A[0]);
			push (@cdt_lines, $_);
	        
            
        }

	# close in file handle
        close(IN);
	closedir(DIR);
        
        # store size of bin array
        $size_A = @cdt_id;


print "Contains: $size_A data rows\n";


################################################################################
# Read in the B file values to four arrays
################################################################################


opendir(DIR,$inB_indir_path) || die "Unable to access file at: $inB_indir_path $!\n";

@files_B = readdir(DIR);

# process the input file within indir_path
foreach $infile_B (@files_B){    

    # ignore hidden files and only get those ending .txt
    if (($infile_B !~ /^\.+/) && ($infile_B =~ /.*\.txt/)){
        
       

print "Found, and processing, $infile_B \n";

open(IN, "$inB_indir_path/$infile_B")
            || die "Unable to open $infile_B: $!";
        
        # define five new arrays to store required values from infile
        my @C3_id; # the C3.txt file ORF ID string
        my @C3_lines; # the c3.txt file full data lines
      
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_B = split('\t',$_);

            # store the columns we want in the new arrays
            push(@C3_id,$line_B[0]);
			push(@C3_lines,$_);
            
        }
        
        # close in file handle
        close(IN);
	closedir(DIR);

# store size of bin array
        $size_B = @C3_id;

print "Contains: $size_B data rows\n";

# check that the .cdt file has one more line than the c3.txt file
if ($size_B != $size_A-1){

	print "Oooh err, your .cdt and c3.txt files don't seem to have matching rows\n";
	exit;

}


######################################################################################
# The output file
######################################################################################

# define outfile name and set ending to .cdt
        $outfile = $infile_B;
        $outfile .= '_converted.cdt';

# try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        


# some counter variables
my $cdt_count = 2; 
my $C3_count = 1;

my @C3_split; # an array to hold the individual C3 row values
my $C3_cols; # the number of C3 columns
my $i=1; # an incrementor variable


# print the first two line of the .cdt output file 

print (OUT $cdt_lines[0]."\n".$cdt_lines[1]."\n");

# spool through the .cdt list and list ID matches with the C3.txt list
until ($cdt_count == $size_A){ #until 1

	until($C3_count == $size_B){ #until 2

		# search for matches
		if ($cdt_id[$cdt_count] eq $C3_id[$C3_count]) { #if 1
		
		# split C3.txt line by tab delimiter to release the individual data values
		@C3_split = split('\t',$C3_lines[$C3_count]);
		$C3_cols = @C3_split;

		
		
		# output a match line to the .cdt file
		
		# print first 3 columns of .cdt row
		print (OUT $cdt_id[$cdt_count]."\t".$cdt_id[$cdt_count]."\t"."1");
		
		#print the C3 data values
		
		for ($i=1; $i <= $C3_cols-1; $i++){
		print(OUT "\t".$C3_split[$i]);
		}
		
		# end the line
		print (OUT "\n");
		
		$C3_count ++;
		

		
		} #if 1  closer
		else{
		$C3_count ++;

		}
	
	}#until 2 closer
		
		$C3_count  =1;
		$cdt_count ++;

	
} #until 1 closer


print "Have just created $outfile\n";
        # close out file handle
        close(OUT);

}}}}
       
