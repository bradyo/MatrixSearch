#!/usr/bin/perl -w
use constant USAGE =><<END;
SYNOPSIS:
fastaToMatrix.pl <fasta input> <html output>

DESCRIPTION:
Converts a multiple sequence fasta file to an html file defining a search
matrix that can be used by the SearchMatrix.pm module. This program takes the
sequences in the Fasta file, aligns them using ClustalW, calculates a
frequency weight vector (ci), and stores the alignment and frequency matrix
in the given html output file.

OPTIONS:
<fasta input>
    input file in fasta format with multiple sequences for alignemnt
<html output>
    outputs the resulting matrix to the given file

EXAMPLES:

fastaToMatrix.pl ctcfSeqs.fa ctcf.html
 
KNOWN ISSUES:
Requires bioperl to be installed.
Requires the environment variable CLUSTALDIR to be set to the location of
the clustalw executable. 
 
AUTHOR:
Brady Olsen
END

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::SimpleAlign;

# get parameters
die(USAGE) if (@ARGV < 2);
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

# load in sequence objects from input fasta file
my $seqIn = Bio::SeqIO->new(-file => "<$inputFilename", -format => 'Fasta');
my @seqObjects = ();
while (my $seq = $seqIn->next_seq()) {
	push(@seqObjects, $seq);
}

# align the sequences using ClustalW
my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM', 'quiet' => 1);
my $clustal = Bio::Tools::Run::Alignment::Clustalw->new(@params);
my $simpleAlign = $clustal->align(\@seqObjects);
die ("Alignement failed\n") if (!defined($simpleAlign));
my @seqs = ();
foreach my $seq ($simpleAlign->each_seq()) {
	my $seqStr = $seq->seq();
	$seqStr =~ s/\./-/g; # replace '.' with '-'
	push(@seqs, $seqStr);
}
my $numSequences = @seqs;
my $numPositions = length($seqs[0]);

# create frequency matrix
my @freqMatrix = ();
for my $pos (1 .. $numPositions) {
	foreach my $base ('A', 'C', 'G', 'T') {
		$freqMatrix[$pos]{$base} = 0;
	}
	for my $seq (@seqs) {
		my $base = substr($seq, $pos - 1, 1);
		$freqMatrix[$pos]{$base}++;	
	}
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
print OUTFILE "\nSequence Alignment:<br>";
print OUTFILE "\n<table name='alignmentTable' border='1'>";
print OUTFILE "\n<tr><td><b>Id</b><td><b>Sequence</b>";
foreach my $seq ($simpleAlign->each_seq()) {
	print OUTFILE "\n<tr><td>", $seq->id, "<td>", $seq->seq;
}
print OUTFILE "\n</table>";
print OUTFILE "\n<br><br>";
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
close OUTFILE;