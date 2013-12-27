#!/usr/bin/perl -w

######################################################################
#
# shuffle-aln.pl --- Randomization of a multiple sequence alignment
#
# Stefan Washietl <wash@tbi.univie.ac.at>
#
# Time-stamp: <04/03/10 13:32:45 wash>
#
######################################################################

use strict "vars";
use Getopt::Long;

######################################################################
#
# void usage();
#
# Print manual page and exit.
#
######################################################################

sub usage{
print <<EOF;

NAME

  shuffle-aln.pl - Randomization of a multiple sequence alignment

DESCRIPTION

  A multiple sequence alignment is read from STDIN and a shuffled
  version is written to STDOUT. FASTA and ClustalW formats are
  supported.

USAGE

   shuffle-aln.pl [OPTIONS] < alignment.aln

OPTIONS

  -f, --format Output format [clustal|fasta]. Default is input format.

  -m, --mode   Mode of shuffling [conservative|complete]. If set to
               "conservative" only columns with the same gap-pattern
               and conservation pattern are shuffled. "complete"
               shuffles all columns. Default is "conservative".

  -n, --n      If a integer n>1 is given here, n shuffled versions
               are written to files named "sampleX.aln"
               (see --mask). Default is 1 and output to STDOUT.

  --mask       Template for output file names. Needs one X
               which is replaced by the current sample number. Default
               is "sampleX.aln".

  -h, --help   Display this help message.
EOF
exit;
}


######################################################################
#
# @shuffledAln shuffleAlnConservative(\@aln)
#
# Shuffles randomly columns of the alignment that have the same
# gap pattern and local conservation pattern.
#
# Takes reference to a array of array (rows/columns) which holds
# the alignment and returns the shuffled array.
#
######################################################################

sub shuffleAlnConservative{

  my @aln=@{$_[0]};
  my $maxRow=$#aln;
  my $maxCol=@{$aln[0]}-1;

  my @list=(); # stores mask for each column
  my %hash=(); # stores columns for each mask as hash of array
               # with mask as key

  # creates characteristic mask for each column and
  # writes @list and %hash
  foreach my $currCol (0..$maxCol){
	my %seen=();
	my $mask='';
	my $counter=0;
	foreach my $currRow (0..$maxRow){
	  my $currNt=$aln[$currRow][$currCol];
	  if ($currNt eq '-'){
		$mask.='-';
	  } else {
		if ($seen{$currNt}){
		  $mask.=$seen{$currNt};
		} else {
		  $counter++;
		  $seen{$currNt}=$counter;
		  $mask.=$counter;
		}
	  }
	}
	push @list, $mask;
	if (!exists $hash{$mask}){
	  $hash{$mask}=[$currCol];
	} else {
	  push @{$hash{$mask}},$currCol;
	}
  }

  # each list of columns with the same mask
  # are shuffled (Fisher-Yates, code from perlfaq4)
  foreach my $arrayRef (values %hash){
	my $i;
	for ($i = @$arrayRef; --$i; ) {
	  my $j = int rand ($i+1);
	  @$arrayRef[$i,$j] = @$arrayRef[$j,$i];
	}
  }

  # columns are reassembled to a shuffled alignment
  my @shuffledAln;
  foreach my $currCol (0..$maxCol){
	my $randomCol=shift @{$hash{$list[$currCol]}};
	foreach  my $currRow (0..$maxRow){
	  $shuffledAln[$currRow][$currCol]=$aln[$currRow][$randomCol];
	}
  }

  return @shuffledAln;
}

######################################################################
#
# @shuffledAln shuffleAlnComplete(\@aln)
#
# Shuffles randomly _all_ columns of the alignment.
#
# Takes reference to a array of array (rows/columns) which holds
# the alignment and returns the shuffled array.
#
######################################################################


sub shuffleAlnComplete{
  my @aln=@{$_[0]};
  my $maxRow=$#aln;
  my $maxCol=@{$aln[0]}-1;
  my $i;
  for ($i = $maxCol+1; --$i; ) {
	my $j = int rand ($i+1);
	for my $k (0..$maxRow){
	  my $tmp=$aln[$k][$j];
	  $aln[$k][$j]=$aln[$k][$i];
	  $aln[$k][$i]=$tmp;
	}
  }
  return @aln;
}

######################################################################
#
# void writeAln(\@aln, $format, [$filename])
#
# Writes alignment given as reference to array of array to the file
# given by $filename. If no $filename is specified writes to STDOUT.
# The alignment format is specified by $format ('clustal' or 'fasta').
#
######################################################################

sub writeAln{
  my @aln=@{$_[0]};
  my $format=$_[1];
  my $file=$_[2];
  my $fileHandle;
  if ($file){
	open(OUTFILE,">$file");
	$fileHandle=OUTFILE;
  } else {
	$fileHandle=STDOUT;
  }

  my $maxRow=$#aln;
  my $maxCol=@{$aln[0]}-1;

  if ($format eq 'fasta'){
	for my $i (0..$maxRow){
	  print $fileHandle ">seq$i\n";
	  my $string='';
	  for my $j (0..$maxCol){
		$string.=$aln[$i][$j];
	  }
	  $string =~ s/(.{60})/$1\n/g;
	  print $fileHandle $string;
	  print $fileHandle "\n";
	}
	print $fileHandle "\n";
	close(OUTFILE) if ($file);
  }

  if ($format eq 'clustal'){
	my $maxName=0;
	my @names=();
	push @names, "seq$_" foreach (0..$maxRow);
	foreach my $name (@names){
	  $maxName=($maxName<length($name))?length($name):$maxName;
	}

	for my $i (0..$#names){
	  my $buffer=" "x(($maxName+6)-length($names[$i]));
	  $names[$i].=$buffer;
	}
	#my $columnWidth=60;
	my $columnWidth=100000000;
	my $currPos=0;

	my @seqs=();
	push @seqs, join('',@$_) foreach (@aln);
	
	my $length=length($seqs[0]);
	print $fileHandle "CLUSTAL W(1.81) multiple sequence alignment\n\n\n";
	while ($currPos<$length){

	  for my $i (0..$#names){
		
		print $fileHandle $names[$i];
		print $fileHandle substr($seqs[$i],$currPos,$columnWidth);
		print $fileHandle "\n";
	  }
	  print $fileHandle "\n\n";
	  $currPos+=$columnWidth;
	}
  }
}

######################################################################
# MAIN PROGRAM
######################################################################

# Get user options

my $outFormat='fasta'; # will be set to detected input format if no --format is given
my $format='';
my $mode='conservative';
my $mask='sampleX.aln';
my $n=1;
my $help='';

# print help message also for '-help' and '?' (not covered by Getopt)
foreach (@ARGV){
  if (($_ eq '-help') or ($_ eq '?')){
	usage();
  }
}

GetOptions ("format:s" => \$format,
			"f:s" => \$format,
			"mask:s" => \$mask,
			"n:s"=>\$n,
			"mode:s" => \$mode,
			"m:s" => \$mode,
			"help"=>\$help,
			"h"=>\$help
		   );

usage() if ($help);

$mode=lc($mode);
$format=lc($format);

if (($mode ne 'complete') and ($mode ne 'conservative')){
  print "Unknown shuffling mode. Use \"complete\" or \"conservative\".\n";
  exit;
}

if (($format ne 'clustal') and ($format ne 'fasta') and ($format)){
  print "Unknown output format. Use \"clustal\" or \"fasta\".\n";
  exit;
}


# Read input from STDIN
my @aln=();
my $currSeq;
while (my $line=<>){
  next if ($line=~/^\s+$/);
  # starting with '> seqname", obviously fasta formatted
  if ($line=~/\s*>\s*(.*$)/){
	$outFormat='fasta';
	while ($line=<>){
	  if ($line=~/\s*>\s*(.*$)/){
		push @aln,[split(//,$currSeq)] if ($currSeq);
		$currSeq="";
		next;
	  }
	  chomp($line);
	  $currSeq.=$line;
	}
	push @aln,[split(//,$currSeq)];
	last;
  }

  # a line with "Clustal" indicates a clustal format, code partly
  # from Bioperl
  my %order;
  my $order;
  my %alignments;
  if ($line=~/clustal/i){
	$outFormat='clustal';
	while(<>) {
	  next if ( /^\s+$/ );	
	  my ($seqname, $aln_line) = ('', '');	
	  if( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ ) {
		($seqname,$aln_line) = ("$1/$2-$3",$4);
	  } elsif( /^(\S+)\s+([A-Z\-]+)\s*$/ ) {
		($seqname,$aln_line) = ($1,$2);
	  } else {
		next;
	  }

	  if( !exists $order{$seqname} ) {
		$order{$seqname} = $order++;
	  }
	  $alignments{$seqname}.= $aln_line;
	}
	foreach my $name ( sort {$order{$a} <=> $order{$b}} keys %alignments){
	  push @aln, [split(//,$alignments{$name})];
	}
  }
}

$outFormat=$format if ($format);

# write a single shuffled version to STDOUT
if ($n==1){
  my @shuffledAln=();

  if ($mode eq 'conservative'){
	@shuffledAln=shuffleAlnConservative(\@aln);
  } else {
	@shuffledAln=shuffleAlnComplete(\@aln);
  }
  writeAln(\@shuffledAln,$outFormat);

# write n shuffled versions to files
} elsif ($n>1){
  for my $i (1..$n){
	my $fileName=$mask;
	if (!($fileName=~s/X/$i/)){
	  die("Invalid mask format (eg. sampleX.aln)");
	}
	my @shuffledAln=();
	if ($mode eq 'conservative'){
	  @shuffledAln=shuffleAlnConservative(\@aln);
	} else {
	  @shuffledAln=shuffleAlnComplete(\@aln);
	}
	writeAln(\@shuffledAln,$outFormat,$fileName);
  }
}
