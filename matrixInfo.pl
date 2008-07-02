use constant USAGE =><<END;
SYNOPSIS:
matrixInfo.pl <matrix html file>

DESCRIPTION:
This program takes in an html file describing a frequency matrix generated
by the matrix creation program (align_convert.pl). This program requires the
SearchMatrix.pm perl module to run.

EXAMPLES:
# load and display the matrix data in matrix.html
matrixInfo.pl matrix.html
 
AUTHOR:
Brady Olsen
END

use SearchMatrix;
use strict;

my $matrixFilename = $ARGV[0];

# check parameters
if (!defined($ARGV[0])) {
	die(USAGE);
}
if (!-e $matrixFilename) {
	die("File not found: $matrixFilename\n");
}

my $searchMatrix = SearchMatrix->new($matrixFilename, $matrixFilename);
$searchMatrix->print();
