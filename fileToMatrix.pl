#!/usr/bin/perl -w
use constant USAGE =><<END;
SYNOPSIS:
fileToMatrix.pl <data input> <html output>

DESCRIPTION:
Converts a data file containing base frequencies into an html file defining
a search matrix that can be used by the SearchMatrix.pm module. This program
takes the frequency data, calculates a frequency weight vector (ci), and stores
frequency matrix in the given html output file.

Input data file should have a 1 line header, followed by any number of rows
(one for each position in the matrix). Each row should contain 4 numbers,
representing relative frequencies of A, C, G, and T respectively. Example:
A       C       G       T        
0.5     0.4     0.0     0.1
0.2     0.6     0.2     0.0

OPTIONS:
<fasta input>
    input file in fasta format with multiple sequences for alignemnt
<html output>
    outputs the resulting matrix to the given file

EXAMPLES:
fileToMatrix.pl ctcf.dat ctcf.html
 
AUTHOR:
Brady Olsen
END

use strict;

# get parameters
die(USAGE) if (length(@ARGV) < 2);
my $inputFilename = $ARGV[0];
my $outputFilename = $ARGV[1];

# check parameters
if (!defined($inputFilename)) {
	die("File not found: $inputFilename\n");
}
if (-e $outputFilename) {
	# warn user if output file already exists
	print("Output file already exists: $outputFilename\n");
	print("Would you like to overwrite? (y\\n): ");
	chomp(my $ans = <STDIN>);
	die("Canceled\n") if ($ans !~ m/^\s*y\s*$/i) # not 'y' or 'Y'
}

# read data file into frequency matrix
open(INFILE, $inputFilename)
    || die("Could not open file: $inputFilename\n");
my @freqMatrix = ();
my $pos = 1;
my $header = <INFILE>; # discard header
while (my $line = <INFILE>) {
	chomp($line);
	my @data = split(/\s+/, $line); # split at white space
	$freqMatrix[$pos]{'A'} = $data[0];
	$freqMatrix[$pos]{'C'} = $data[1];
	$freqMatrix[$pos]{'G'} = $data[2];
	$freqMatrix[$pos]{'T'} = $data[3];
	$pos++;
}
my $numPositions = $pos - 1;
close(INFILE);

# compute Ci array as defined by the Genomatix algorithm
my @ci = ();
for my $pos (1 .. $numPositions) {
	my $sum = 0;
	foreach my $base ('A', 'C', 'G', 'T') {
		my $relativeFreq = $freqMatrix[$pos]{$base};
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