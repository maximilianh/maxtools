#!/usr/bin/perl 
##############################################################################
# Syntax  : perl fastaFilter.pl -i input_Fasta_file -e file_with_sequence_Defs_to_be_excluded -o filtered.fas 
#
# Note: This script is used to filter out sequences.
#
# Written by Philippe Chouvarine, Mississippi Genome Exploration Laboratory.
#
###############################################################################


use Getopt::Std;

getopt ('ioe');

open (FD_E, "$opt_e") || die "Error opening $opt_e";
open (FD_IN, "$opt_i") || die "Error opening $opt_i";
open(FD_OUT, ">$opt_o") || die("could not write newseq file");


%exclude = (); # this hash will contain records to be exluded for fast comparison

while($eline = <FD_E>)
{
	$exclude{$eline} = 1;
}
close(FD_E);

while($line = <FD_IN>) 
{
	if($line =~ /^(\s*>(\S+).*)$/)   # >line
	{
		$badseq = 0;			

		$content = substr($line, 1, length($line));  # >line

		if($exclude{$content} == 1) 
		{
			$badseq = 1;
			print "$content\n";
		}
	}
			
	
	if($badseq==0)
	{
		print FD_OUT $line;
	}
	

	
}

close(FD_OUT);
close(FD_IN);

