#!/usr/bin/perl -w

=head1 NAME

matrixSearch - program to do a similarity matrix search on a fasta file

=head1 SYNOPSIS

matrixSearch.pl <input fasta> <output> [OPTIONS]

=head1 DESCRIPTION

Searches a given fasta file for matches using matrix definitions generated
by the fastaToMatrix.pl program. The program searches using the Genomatix
MatInspector algorithm, using core similairty and matrix similairty scores.
Requres SearchMatrix.pm.

=head1 OPTIONS

<input fasta>
    fasta file containing the sequence or sequences to search
<output>
    file to write the comma separated data to
-matrix = <string>
    the html matrix files (generated by align_convert.pl) to seach input for.
    Multiple matricies must be separated by a comma with no spaces.
-coreSim = <float> (Default = 0)
    Only output results exceeding the specified core similarity threshold.
	Must be between 0 and 1.
-matSim = <float> (Default = 0.8)
    Only output results exceeding the specified matrix similarity threshold.
	Must be between 0 and 1.

This section intentionally left blank.

=head1 EXAMPLES

matrixSearch.pl chr7.fa output.txt
    -matrix=lm7.html,lm2.html,lm22.html
    -coreSim=1 -matSim=0.8

=head1 AUTHOR

Brady Olsen

=cut

use Getopt::Std;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Benchmark;
use SearchMatrix;
use strict;

use constant USAGE => "usage";

# get command line parameters
die(USAGE) if (@ARGV < 2);
my $inputFilename = $ARGV[0];
my $outputFilename = $ARGV[1];
my $matrix = "";
my $coreSimilarity = 0;
my $matrixSimilarity = 0.8;
GetOptions(
    "matrix=s" => \$matrix,
    "coreSim=s" => \$coreSimilarity,
    "matSim=s" => \$matrixSimilarity)
    or die($!);
my @matrixFilenames = split(',', $matrix);

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
if (@matrixFilenames < 1) {
    die("No search matrix specified\n");
} else {
    foreach my $matrixFilename (@matrixFilenames) {
        if (!-e $matrixFilename) {
            die("Matrix file not found: $matrixFilename");
        }
    }
}
if ($coreSimilarity < 0 || $coreSimilarity > 1) {
    die("coreSim must be between 0 and 1\n");
}
if ($matrixSimilarity < 0 || $matrixSimilarity > 1) {
    die("matSim must be between 0 and 1\n");
}

# load matricies
my @searchMatricies = ();
foreach my $matrixFilename (@matrixFilenames) {
    my $matrix = new SearchMatrix($matrixFilename, $matrixFilename);
    push(@searchMatricies, $matrix);
}

# get max matrix length from search matricies
my $maxMatrixLength = 0;
foreach my $matrix (@searchMatricies) {
    if ($maxMatrixLength < $matrix->{_matrixLength}) {
        $maxMatrixLength = $matrix->{_matrixLength};
    }
}

# load sequences from command line parameter
my $in = Bio::SeqIO->new(-file => "<$inputFilename", -format => 'Fasta')
    || die("Could not load $inputFilename with Bioperl: $!\n");
my @seqObjects = ();
while (my $seq = $in->next_seq()) {
	push(@seqObjects, $seq);
}

# scan each input sequence with each matrix and output to file
open(FILEOUT, ">$outputFilename")
    || die("Could not open file: $outputFilename\n");
print FILEOUT "matrix, sequence, start, end, strand, score, sequence\n";
    
my $startTime = new Benchmark();
foreach my $seqObject (@seqObjects) {
	my $sequenceLength = length($seqObject->seq);
	
	$| = 1; # flush output after each print
	print "Seq: ", $seqObject->id, "\n";
	
	my $buffer = $seqObject->seq;
	for (my $pos = 1; $pos <= $sequenceLength; $pos++) {	
		foreach my $searchMatrix (@searchMatricies) {
			my $matrixLength = $searchMatrix->{_matrixLength};
			if ($pos + $matrixLength - 1 <= $sequenceLength) {
				# process forward strand
				my $subSeq = substr($buffer, 0, $matrixLength);
				my $coreSim = $searchMatrix->getCoreSimilarity($subSeq);
				if ($coreSim >= $coreSimilarity) {
					my $matrixSim = $searchMatrix->getSimilarity($subSeq);			
					if ($matrixSim >= $matrixSimilarity) {
						# output the match
						my $end = $pos + $matrixLength;
						print FILEOUT "$searchMatrix->{_name}, ",
                            $seqObject->id, ", $pos, $end, +, ",
							"$matrixSim, $subSeq\n"; 
					}
				}
	
				# process reverse complement strand
				$subSeq = reverse($subSeq);
				$subSeq =~ tr/ACGT/TGCA/;
				$coreSim = $searchMatrix->getCoreSimilarity($subSeq);
				if ($coreSim >= $coreSimilarity) {
					my $matrixSim = $searchMatrix->getSimilarity($subSeq);
					if ($matrixSim >= $matrixSimilarity) {
						# output the match
						my $end = $sequenceLength - $pos + 1;
						my $start = $end - $matrixLength + 1;
						print FILEOUT "$searchMatrix->{_name}, ",
                            $seqObject->id, ", $start, $end, -, ",
							"$matrixSim, $subSeq\n"; 
					}
				}
			}
		}
	
		# shrink buffer
		$buffer = substr($buffer, 1);
	}
}
close(FILEOUT);

my $endTime = new Benchmark;
my $deltaTime = timediff($endTime, $startTime);
print "Total Time = ", timestr($deltaTime), "\n";

