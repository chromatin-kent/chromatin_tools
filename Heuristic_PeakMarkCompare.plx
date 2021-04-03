#!/usr/bin/perl
# Written: Nick Kent, 14th June 2011
# Last updated: Nick Kent, 27th Nov 2013
# Last updated: Nick Kent, 5th Sep 2019 to use Cwd
# USAGE:- perl PeakMarkCompare.plx
#
#
# This script takes two .sgrs (termed A and B) representing a particular Partn value and uses various
# criteria to mark peaks in both and to describe those peak bins that match and do not match
# between the data sets.
#
#
# NOTE: The input files for this script must be ORDERED IDENTICALLY; so the 
# number of bins per chromosome and all chromsosomes are in the same order. The script
# does a simple block-by-block pair-wise comparison between the files in order to run
# quickly. There is little error checking at the moment for files that do not match up, and
# the script will therefore vomit complete nonsense if you get this wrong.
#
# The script begins with a simple peak marking process which calls peak summit bins
# and sets all other bins to zero. The summit dyad read frequency of the B input .sgr file
# can be scaled by setting the $scale_factor variable. The number here should be based on some
# estimate of the difference in read-depth between the A file and B file data. The most
# sensible derivation of this would seem to be to use the CFD sum values from a SiteWriter.pl
# analysis for this particular Partn or the actual experimental read depth.
#
# In order to filter out low level noise, the script applies a read theshold  which peak
# summits must exceed (set with the $read_thresh variable). However, to avoid missing peaks
# which lie just at the threshold, the script also gives a boundary score to bins which
# contain peaks occuring within a set percentage of the read threshold regardless of whether 
# or not they get marked. When the script checks for peaks which do NOT match between the
# A and B arrays it then also checks for the presence of a boundary condition in either data
# set at that point. If one occurs, then any non-matching peak is ignored and not sent to the 
# NOT output file. This stops a peak in one data set from appearing as unique just because
# its partner in the other file was just below the read threshold. Applying this correction is
# a crude but effective way of cutting down erroneous peak calls.
#
# In addition to making yes/no peak calls, this script also uses the scaled read 
# frequency information to decide further NOT conditions based on a user-defined read frequency 
# difference. i. e. the fold freq difference is greater than 2 or whatever (set with
# the $diff variable). NOTE in this version of the script such read frequency differences
# are noted in the NOT output file and are deleted from the MATCHes - this is different
# to the other versions of this script
#
# The script outputs 2 .sgr files showing A MATCH B and A NOTMATCH B summit bins.
# Note: The MATCH logic is relatively simple - all the clever thresholding is applied
# to the NOT files (which attempt to show the interesting differences!). The MATCH file
# therefore probably represents an understimate.
#
# These files can then be further processed using Interval scoring scripts to select for those
# occuring in regions of interest, such as ORFs or intergenic regions.
#
# WARNING: this is still a development script;
################################################################################

use strict;
use warnings;
use Math::Round;
use Cwd;

################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $sgr_A_indir_path   - The directory containing the full genome Partn.sgr file A
# $sgr_B_indir_path   - The directory containing the full genome Partn.sgr file B
# $outdir_path  - The directory to store the .sgr and .cfd output files
# $bin_window - 5 or 15bp (should have called it $bin_size, but what the hell.
# $diff -  fold difference between read frequencies to trigger a NOT peak.
# $read_thresh  - The aligned read number noise threshold value
# $scale_factor - A proportion based on differences in read depth to be applied to
# the B file.
# $bound_window - The percentage above and below $read_thresh at which to mark a
# peak or bin as a "boundary". Peaks in this category will NOT be counted as 
# differences between data sets (0.1 = 10% - suggest this is set to 0.3)
################################################################################

my $sgr_A_indir_path =cwd."/A_in";
my $sgr_B_indir_path =cwd."/B_in";
my $outdir_path =cwd."/PMC_files_out";
my $descriptor;
my $bin_window = 10; 
my $diff = 1.5;
my $read_thresh = 15;
my $scale_factor = 1.04;
my $bound_window = 0.3;


################################################################################
################################################################################
#####################          MAIN PROGRAM          ###########################
################################################################################
################################################################################

# define some variables

my $infile_sgr_A;
my $infile_sgr_B;
my $outfile;
my @line_sgr_A;
my @line_sgr_B;
my @files_sgr_A;
my @files_sgr_B;
my $sgr_A_size;
my $sgr_B_size;
my $A_counter=1;
my $B_counter=1;
my $counter=0;
my $out_freq =10;# an arbitrary value for the output .sgr file frequencies
my $count = 0;
my $no = 0;
my $up_bound = $read_thresh*(1+$bound_window);
my $down_bound = $read_thresh*(1-$bound_window);

################################################################################
# Read in the .sgr A file values to three arrays
################################################################################
print "start:\t", `date`."\n";
print "Read noise threshold is set to: $read_thresh\n";
print "B file scale factor is set to: $scale_factor\n";
print "The threshold boundary is set at: $bound_window times either side of the noise threshold\n";
print "Peak fold-difference trigger is set to: $diff\n\n";

opendir(DIR,$sgr_A_indir_path) || die "Unable to access file at: $sgr_A_indir_path $!\n";

@files_sgr_A = readdir(DIR);

# process the input file within sgr_A_indir_path
foreach $infile_sgr_A (@files_sgr_A){   

# deleted cfd arrays here - might need to put back 

    # ignore hidden files and only get those ending .sgr
    if (($infile_sgr_A !~ /^\.+/) && ($infile_sgr_A =~ /.*\.sgr/)){
        
       

print "Found, and processing, $infile_sgr_A \n";

open(IN, "$sgr_A_indir_path/$infile_sgr_A")
            || die "Unable to open $infile_sgr_A: $!";
        
        # define arrays to store the .sgr values from infile and a threshold border marker
        my @sgr_A_chr;
        my @sgr_A_bin;
	my @sgr_A_freq_in; # the actual freq values
	my @sgr_A_freq; # the remaining freq values after peak marking
	my @sgr_A_threshborder;

        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_sgr_A = split('\t',$_);

            # store the columns we want in the new arrays
            push(@sgr_A_chr,$line_sgr_A[0]);
            push(@sgr_A_bin,$line_sgr_A[1]);
            push(@sgr_A_freq_in,$line_sgr_A[2]);
            push(@sgr_A_freq,$line_sgr_A[2]);
            push(@sgr_A_threshborder,0);
        }
        
        # close in file handle
        close(IN);
		closedir(DIR);

# store size of bin array
        $sgr_A_size = @sgr_A_chr;

print "which contains: $sgr_A_size bin values\n\n";

#
################################################################################
# Read in the .sgr B file values to three arrays
################################################################################


opendir(DIR,$sgr_B_indir_path) || die "Unable to access file at: $sgr_B_indir_path $!\n";

@files_sgr_B = readdir(DIR);

# process the input file within sgr_B_indir_path
foreach $infile_sgr_B (@files_sgr_B){   

# deleted cfd arrays here - might need to put back 

    # ignore hidden files and only get those ending .sgr
    if (($infile_sgr_B !~ /^\.+/) && ($infile_sgr_B =~ /.*\.sgr/)){
        
       

print "Found, and processing, $infile_sgr_B \n";

open(IN, "$sgr_B_indir_path/$infile_sgr_B")
            || die "Unable to open $infile_sgr_B: $!";
        
        # define arrays to store the .sgr values from infile and a threshold border marker
        my @sgr_B_chr;
        my @sgr_B_bin;
	my @sgr_B_freq_in;# the actual freq values
	my @sgr_B_freq;# the remaining freq values after peak marking
	my @sgr_B_threshborder;

        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_sgr_B = split('\t',$_);

            # store the columns we want in the new arrays (APPLIES $scale_factor
	    # to B file read frequency value)
            push(@sgr_B_chr,$line_sgr_B[0]);
            push(@sgr_B_bin,$line_sgr_B[1]);
            push(@sgr_B_freq_in,round($line_sgr_B[2]*$scale_factor));
            push(@sgr_B_freq,round($line_sgr_B[2]*$scale_factor));
            push(@sgr_B_threshborder,0);
        }
        
        # close in file handle
        close(IN);
		closedir(DIR);

# store size of bin array
        $sgr_B_size = @sgr_B_chr;

print "which contains: $sgr_B_size bin values\n\n";




if ($sgr_A_size != $sgr_B_size){
print "Those are not the same size! Suicide is the only option\n";
exit;
}
#########################################################################
# Refine main peaks and those close to the noise threshold in the A file
##########################################################################
        # need a variable to store line count
        $count = 2;
        
        # this calls the peaks - giving an x-axis bin value ONLY for the peak centre and
	# a y-axis value as the peak hight scaled to some value proportionate to relative
	# read depth a relevant pair-wise comparison
	# ALTER HERE to define how many bins either side if a central value are tolerated.
	#
	# ALTER HERE to change >= to > to increase "stringency"? e.g
	
	# ($sgr_A_freq_in[$count-2]>$sgr_A_freq_in[$count-3]) &&
	 #     ($sgr_A_freq_in[$count-1]>$sgr_A_freq_in[$count-2]) &&
	  #    ($sgr_A_freq_in[$count]>=$sgr_A_freq_in[$count-1]) && 
	  #    ($sgr_A_freq_in[$count]>=$sgr_A_freq_in[$count+1]) &&
	   #   ($sgr_A_freq_in[$count+1]>$sgr_A_freq_in[$count+2])&&
	   #   ($sgr_A_freq_in[$count+2]>$sgr_A_freq_in[$count+3])
	   #
	   # ($sgr_A_freq_in[$count-1]>$sgr_A_freq_in[$count-2]) &&
	    #  ($sgr_A_freq_in[$count]>$sgr_A_freq_in[$count-1]) && 
	     # ($sgr_A_freq_in[$count]>=$sgr_A_freq_in[$count+1]) &&
	      #($sgr_A_freq_in[$count+1]>$sgr_A_freq_in[$count+2])


	while ($count < $sgr_A_size-2){
	

	if ( 
	     

	      ($sgr_A_freq_in[$count]>$sgr_A_freq_in[$count-1]) && 
	      ($sgr_A_freq_in[$count]>=$sgr_A_freq_in[$count+1])

	     
	      ){

		# test for clear summit above threshold
		if ($sgr_A_freq_in[$count]>=$up_bound){

		  $count ++;
			}
		# test for marginal summit above threshold
		elsif($sgr_A_freq_in[$count]>=$read_thresh && $sgr_A_freq_in[$count]< $up_bound){

		  $sgr_A_threshborder[$count] = 1;
		  $count ++;
			}
		# test for marginal summit below threshold
		elsif($sgr_A_freq_in[$count] < $read_thresh && $sgr_A_freq_in[$count] > $down_bound){

		  $sgr_A_freq[$count] =$no;
		  $sgr_A_threshborder[$count] = 1;
		  $count ++;
			}
		  }

		
		  $sgr_A_freq[$count] =$no;
		  $count ++;

}


##########################################################################
# Refine main peaks and those close to the noise threshold in the B file
##########################################################################
        # need a variable to store line count
        $count = 2;
        
        # this calls the peaks - giving an x-axis bin value ONLY for the peak centre and
	# a y-axis value as the peak hight scaled to some value proportionate to relative
	# read depth a relevant pair-wise comparison
	# ALTER HERE to define how many bins either side if a central value are tolerated.
	#
	# ALTER HERE to change >= to > to increase "stringency"?
	#
	#($sgr_B_freq_in[$count-1]>$sgr_B_freq_in[$count-2]) &&
	  #    ($sgr_B_freq_in[$count]>$sgr_B_freq_in[$count-1]) && 
	   #   ($sgr_B_freq_in[$count]>=$sgr_B_freq_in[$count+1]) &&
	    #  ($sgr_B_freq_in[$count+1]>$sgr_B_freq_in[$count+2])
	      
	#
	


	while ($count < $sgr_B_size-2){

	if ( 
	      ($sgr_B_freq_in[$count]>$sgr_B_freq_in[$count-1]) && 
	      ($sgr_B_freq_in[$count]>=$sgr_B_freq_in[$count+1]) 
	      
	      ){

		# test for clear summit above threshold
		if ($sgr_B_freq_in[$count]>=$up_bound){
		
		  $count ++;
			}
		# test for marginal summit above threshold
		elsif($sgr_B_freq_in[$count]>=$read_thresh && $sgr_B_freq_in[$count]< $up_bound){

		  $sgr_B_threshborder[$count] = 1;
		  $count ++;
			}
		# test for marginal summit below threshold
		elsif($sgr_B_freq_in[$count] < $read_thresh && $sgr_B_freq_in[$count] > $down_bound){

		  $sgr_B_freq[$count] =$no;
		  $sgr_B_threshborder[$count] = 1;
		  $count ++;
			}
		  }

		
		  $sgr_B_freq[$count] =$no;
		  $count ++;
	

}




#########################################################################
# Count through the A file and look for matches and A NOT B mismatches
#########################################################################

# define some arrays:
my @A_match_chr;
my @A_match_bin;
my @B_match_chr;
my @B_match_bin;
my @AnotB_chr;
my @AnotB_bin;
my @BnotA_chr;
my @BnotA_bin;
my $A_match_size;
my $B_match_size;
my $AnotB_size;
my $BnotA_size;
my $diff_counter;



until ($A_counter == $sgr_A_size){ #until 1

      


	  if ($sgr_A_freq[$A_counter] > 0 && $sgr_B_freq[$A_counter-1] > 0){

	      push(@A_match_chr,$sgr_A_chr[$A_counter]);
	      push(@A_match_bin,$sgr_A_bin[$A_counter]);

		  # then test to see if A/B dyad read frequency is greater than $diff
		  # if it then this isn't a real match so remove ot from the MATCH list
		  if (($sgr_A_freq[$A_counter]/$sgr_B_freq[$A_counter-1]) > $diff){

		  push(@AnotB_chr,$sgr_A_chr[$A_counter]);
		  push(@AnotB_bin,$sgr_A_bin[$A_counter]);
		  pop @A_match_chr;
		  pop @A_match_bin;
		  $diff_counter ++;
		  }
	     $A_counter ++;

	   }elsif ($sgr_A_freq[$A_counter] > 0 && $sgr_B_freq[$A_counter] > 0){

	      push(@A_match_chr,$sgr_A_chr[$A_counter]);
	      push(@A_match_bin,$sgr_A_bin[$A_counter]);

		  # then test to see if A/B dyad read frequency is greater than $diff
		  # if it then this isn't a real match so remove ot from the MATCH list
		  if (($sgr_A_freq[$A_counter]/$sgr_B_freq[$A_counter]) > $diff){

		  push(@AnotB_chr,$sgr_A_chr[$A_counter]);
		  push(@AnotB_bin,$sgr_A_bin[$A_counter]);
		  pop @A_match_chr;
		  pop @A_match_bin;
		  $diff_counter ++;
		  }
	     $A_counter ++;

	   }elsif ($sgr_A_freq[$A_counter] > 0 && $sgr_B_freq[$A_counter+1] > 0){

	      push(@A_match_chr,$sgr_A_chr[$A_counter]);
	      push(@A_match_bin,$sgr_A_bin[$A_counter]);

		  # then test to see if A/B dyad read frequency is greater than $diff
		  # if it then this isn't a real match so remove ot from the MATCH list
		  if (($sgr_A_freq[$A_counter]/$sgr_B_freq[$A_counter+1]) > $diff){

		  push(@AnotB_chr,$sgr_A_chr[$A_counter]);
		  push(@AnotB_bin,$sgr_A_bin[$A_counter]);
		  pop @A_match_chr;
		  pop @A_match_bin;
		  $diff_counter ++;
		  }
	      $A_counter ++;

	    # Mark an A NOT B for any A peak as long as there is no B boundary condition present
	    }elsif ($sgr_A_freq[$A_counter] > 0 && ($sgr_B_freq[$A_counter-1]+$sgr_B_freq[$A_counter]+$sgr_B_freq[$A_counter+1]) == 0 &&
	    ($sgr_B_threshborder[$A_counter-1]+$sgr_B_threshborder[$A_counter]+$sgr_B_threshborder[$A_counter+1]) == 0){

		  push(@AnotB_chr,$sgr_A_chr[$A_counter]);
		  push(@AnotB_bin,$sgr_A_bin[$A_counter]);
		  
		  $A_counter ++;

	    # Ignore any A peak with a boundary mark if B boundary mark is also present 
	    }elsif ($sgr_A_freq[$A_counter] > 0 && ($sgr_B_freq[$A_counter-1]+$sgr_B_freq[$A_counter]+$sgr_B_freq[$A_counter+1]) == 0 &&
	    ($sgr_A_threshborder[$A_counter-1]+$sgr_A_threshborder[$A_counter]+$sgr_B_threshborder[$A_counter+1]) > 0 &&
	    ($sgr_B_threshborder[$A_counter-1]+$sgr_B_threshborder[$A_counter]+$sgr_B_threshborder[$A_counter+1]) > 0){
		  
		  $A_counter ++;

	      }else{
	      $A_counter ++;	     
	      }
}

 

$A_match_size = @A_match_chr;

print "\n$A_match_size peaks MATCH between A and B, taking into account dyad frequencies\n";
print "$diff_counter matches showed an A/B frequency ratio greater than $diff and were transferred to the A NOT B list.\n\n";

$AnotB_size = @AnotB_chr;

print "The A NOT B file identifies $AnotB_size peaks\n\n";



####################################################################################
# now do it all again to check for B not A (I'm sure there's a smarter way...?)
####################################################################################


$diff_counter = 0;
$A_counter = 1;

until ($A_counter == $sgr_B_size){ #until 1

	      


	  if ($sgr_B_freq[$A_counter] > 0 && $sgr_A_freq[$A_counter-1] > 0){

	      push(@B_match_chr,$sgr_B_chr[$A_counter]);
	      push(@B_match_bin,$sgr_B_bin[$A_counter]);

		  # then test to see if A/B dyad read frequency is greater than $diff
		  if (($sgr_B_freq[$A_counter]/$sgr_A_freq[$A_counter-1]) > $diff){

		  push(@BnotA_chr,$sgr_B_chr[$A_counter]);
		  push(@BnotA_bin,$sgr_B_bin[$A_counter]);
		  pop @B_match_chr;
		  pop @B_match_bin;
		  $diff_counter ++;
		  }
	     $A_counter ++;

	   }elsif ($sgr_B_freq[$A_counter] > 0 && $sgr_A_freq[$A_counter] > 0){

	      push(@B_match_chr,$sgr_B_chr[$A_counter]);
	      push(@B_match_bin,$sgr_B_bin[$A_counter]);

		  # then test to see if A/B dyad read frequency is greater than $diff
		  if (($sgr_B_freq[$A_counter]/$sgr_A_freq[$A_counter]) > $diff){

		  push(@BnotA_chr,$sgr_B_chr[$A_counter]);
		  push(@BnotA_bin,$sgr_B_bin[$A_counter]);
		  pop @B_match_chr;
		  pop @B_match_bin;
		  $diff_counter ++;
		  }
	     $A_counter ++;

	   }elsif ($sgr_B_freq[$A_counter] > 0 && $sgr_A_freq[$A_counter+1] > 0){

	      push(@B_match_chr,$sgr_B_chr[$A_counter]);
	      push(@B_match_bin,$sgr_B_bin[$A_counter]);

		  # then test to see if A/B dyad read frequency is greater than $diff
		  if (($sgr_B_freq[$A_counter]/$sgr_A_freq[$A_counter+1]) > $diff){

		  push(@BnotA_chr,$sgr_B_chr[$A_counter]);
		  push(@BnotA_bin,$sgr_B_bin[$A_counter]);
		  pop @B_match_chr;
		  pop @B_match_bin;
		  $diff_counter ++;
		  }
	      $A_counter ++;

	    # Mark a B NOT A for any B peak as long as there is no A boundary condition present
	    }elsif ($sgr_B_freq[$A_counter] > 0 && ($sgr_A_freq[$A_counter-1]+$sgr_A_freq[$A_counter]+$sgr_A_freq[$A_counter+1]) == 0 &&
	    ($sgr_A_threshborder[$A_counter-1]+$sgr_A_threshborder[$A_counter]+$sgr_A_threshborder[$A_counter+1]) == 0){

		  push(@BnotA_chr,$sgr_B_chr[$A_counter]);
		  push(@BnotA_bin,$sgr_B_bin[$A_counter]);
		  
		  $A_counter ++;

	    # Ignore any B peak with a boundary mark if A boundary mark is also present 
	    }elsif ($sgr_B_freq[$A_counter] > 0 && ($sgr_A_freq[$A_counter-1]+$sgr_A_freq[$A_counter]+$sgr_A_freq[$A_counter+1]) == 0 &&
	    ($sgr_A_threshborder[$A_counter-1]+$sgr_A_threshborder[$A_counter]+$sgr_B_threshborder[$A_counter+1]) > 0 &&
	    ($sgr_B_threshborder[$A_counter-1]+$sgr_B_threshborder[$A_counter]+$sgr_B_threshborder[$A_counter+1]) > 0){
		  
		  $A_counter ++;

	      }else{
	      $A_counter ++;	     
	      }
}

      
$B_match_size = @B_match_chr;

print "\n$B_match_size peaks MATCH between B and A, taking into account dyad frequencies\n";
print "$diff_counter matches showed an B/A frequency ratio greater than $diff and were transferred to the B NOT A list.\n\n";

$BnotA_size = @BnotA_chr;

print "The B NOT A file identifies $BnotA_size peaks\n\n";


#################################################################################
# Write the outfiles - one MATCH file and one NOTMATCH file that is a concatenation
# of the ANOTB and BNOTA sets
###################################################################################

# define MATCH outfile name from infile names
$outfile = substr($infile_sgr_A,0,-4)."_MATCH_".substr($infile_sgr_B,0,-4)."_t$read_thresh"."_sf$scale_factor"."_b$bound_window"."_d$diff";
$outfile .= '.sgr';

 # try and open output file
open(OUT,"> $outdir_path/$outfile")
|| die "Unable to open $outfile: $!";

$counter =0;

until($counter == $A_match_size){

      print(OUT 	$A_match_chr[$counter]."\t".
			$A_match_bin[$counter]."\t".
			$out_freq."\n");

      $counter ++;
 }

# close out file handle
        close(OUT);

# define NOTMATCH outfile name from infile names
$outfile = substr($infile_sgr_A,0,-4)."_NOTMATCH_".substr($infile_sgr_B,0,-4)."_t$read_thresh"."_sf$scale_factor"."_b$bound_window"."_d$diff";
$outfile .= '.sgr';

 # try and open output file
open(OUT,"> $outdir_path/$outfile")
|| die "Unable to open $outfile: $!";

$counter =0;

until($counter == $AnotB_size){

      print(OUT 	$AnotB_chr[$counter]."\t".
			$AnotB_bin[$counter]."\t".
			$out_freq."\n");

      $counter ++;
 }


$counter =0;

until($counter == $BnotA_size){

      print(OUT 	$BnotA_chr[$counter]."\t".
			$BnotA_bin[$counter]."\t".
			$out_freq."\n");

      $counter ++;
 }

# close out file handle
        close(OUT);
print "end:\t", `date`."\n";
}}}}
