#!/usr/bin/perl -w
use constant USAGE =><<END;
SYNOPSIS:
arrayToMatrix.pl <output html>

DESCRIPTION:
Converts frequency arrays (defined in code) into an html matrix file that
can be used by the SearchMatrix.pm module.

OPTIONS:
<output_file>
    outputs the resulting matrix to the given file

EXAMPLES:
arrayToMatrix.pl matrix.html
 
AUTHOR:
Brady Olsen
END

use strict;

# modify the following variables according to your matrix
my @freqA = (13,4,12,16,4,1,39,3,7,51,0,22,2,1,1,10,23,2,8,18);
my @freqC = (13,20,11,3,49,56,1,31,31,1,0,0,1,0,2,41,1,30,22,16);
my @freqG = (19,11,23,32,2,0,8,23,7,2,57,35,35,56,52,1,29,23,8,19);
my @freqT = (12,22,11,6,2,0,9,0,12,3,0,0,19,0,2,5,4,2,19,4);
my $numSequences = 57;
my $numPositions = length(@freqA);

# check parameter
die(USAGE) if (length(@ARGV) < 1);
my $outputFilename = $ARGV[0];
if (-e $outputFilename) {
	# warn user if output file already exists
	print("Output file already exists: $outputFilename\n");
	print("Would you like to overwrite? (y\\n): ");
	chomp(my $ans = <STDIN>);
	die("Canceled\n") if ($ans !~ m/^\s*y\s*$/i) # not 'y' or 'Y'
}

# create frequency matrix
my @freqMatrix = ();
for my $pos (1 .. $numPositions) {
	$freqMatrix[$pos]{'A'} = $freqA[$pos - 1];
	$freqMatrix[$pos]{'C'} = $freqC[$pos - 1];
	$freqMatrix[$pos]{'G'} = $freqG[$pos - 1];
	$freqMatrix[$pos]{'T'} = $freqT[$pos - 1];
}

# compute Ci array as defined by the Genomatix algorithm
my @ci = ();
for my $pos (1 .. $numPositions) {
	my $sum = 0;
	foreach my $base ('A', 'C', 'G', 'T') {
		my $freq = $freqMatrix[$pos]{$base};
		my $relativeFreq = $freq / $numSequences;
		if ($relativeFreq > 0) {
			$sum += $relativeFreq * log($relativeFreq);
		}
	}
	$ci[$pos] = 100 / log(5) * ($sum + log(5));
}

# output to html file in table format
open(OUTFILE, ">$outputFilename")
    || die("Could not open file: $outputFilename\n");
print OUTFILE "\nNucleotide Distribution Matrix:<br>";
print OUTFILE "\n<table name='matrix' border='1'>";
print OUTFILE "\n<tr><td><b>Pos.</b>";
for my $pos (1 .. $numPositions) {
	print OUTFILE "\n<td><b>&nbsp;&nbsp;", $pos, "</b>";
}
foreach my $base ('A', 'C', 'G', 'T') {
	print OUTFILE "\n<tr name='freq$base'><td><b>$base</b>";
	for my $pos (1 .. $numPositions) {
		print OUTFILE "<td>", $freqMatrix[$pos]{$base};
	}
}
print OUTFILE "\n<tr name='ci'>";
print OUTFILE "<td><b>C<sub>i</sub></b></td>";
for my $pos (1 .. $numPositions) {
	my $ciString = sprintf("%0.1f", $ci[$pos]);
	print OUTFILE "<td>$ciString</td>";
}
print OUTFILE "\n</tr>";
print OUTFILE "\n</table>";
close(OUTFILE);