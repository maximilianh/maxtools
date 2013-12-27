#!/usr/bin/env perl
#
#	fixStepToBedGraph.pl - read fixedStep and variableStep wiggle input data,
#		output four column bedGraph format data
#
#	used thusly:
#	zcat fixedStepData.gz | ./fixStepToBedGraph.pl \
#		| gzip > bedGraph.gz

#	$Id: fixStepToBedGraph.pl,v 1.3 2008/03/05 18:42:01 kate Exp $
#   Wed Nov 24 13:38:27 GMT 2010 track and variableStep changes max

use warnings;
use strict;

my $position = 0;
my $chr = "";
my $step = 1;
my $span = 1;

print STDERR "usage: fixStepToBedGraph.pl\n";
print STDERR "\trun in a pipeline like this:\n";
print STDERR "usage: ",
"zcat fixedStepData.gz | fixStepToBedGraph.pl | gzip > bedGraph.gz",
	"\n";
print STDERR "reading input data from stdin ...\n";

while (my $dataValue = <>)
{
    chomp $dataValue;
    if ( $dataValue =~ m/^track /) {
        print STDERR "Skipping track line\n";
        next;
    }
    if ( $dataValue =~ m/^fixedStep / || $dataValue =~ m/^variableStep / ) {
        $position = 0;
        $chr = "";
        $step = 1;
        $span = 1;
        my @a = split('\s', $dataValue);
        for (my $i = 1; $i < scalar(@a); ++$i) {
            my ($ident, $value) = split('=',$a[$i]);
            if ($ident =~ m/chrom/) { $chr = $value; }
            elsif ($ident =~ m/start/) { $position = $value-1; }
            elsif ($ident =~ m/step/) { $step = $value; }
            elsif ($ident =~ m/span/) { $span = $value; }
            else {
            print STDERR "ERROR: unrecognized fixedStep line: $dataValue\n";
            exit 255;
            }
        }
      } else {
          my @b = split('\s', $dataValue);
          if (scalar(@b)>1) {
              $position = $b[0];
              $dataValue = $b[1];
          }
          printf "%s\t%d\t%d\t%g\n", $chr, $position, $position+$span, $dataValue;
          $position = $position + $step;
      }
}
