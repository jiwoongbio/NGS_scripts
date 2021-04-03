#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	's=s' => \(my $stranded = ''),
);
my ($fastaFile, $motif) = @ARGV;
my ($sequenceName, $sequence) = ('', '');
my $count = 0;
open(my $reader, $fastaFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^>(\S*)/) {
		if($sequence ne '') {
			$count += () = ($sequence =~ /(?=$motif)/g) if($stranded eq '' || $stranded eq 'f' || $stranded eq 'forward');
			$count += () = (getReverseComplementarySequence($sequence) =~ /(?=$motif)/g) if($stranded eq '' || $stranded eq 'r' || $stranded eq 'reverse');
		}
		($sequenceName, $sequence) = ($1, '');
	} else {
		$sequence .= $line;
	}
}
close($reader);
if($sequence ne '') {
	$count += () = ($sequence =~ /(?=$motif)/g) if($stranded eq '' || $stranded eq 'f' || $stranded eq 'forward');
	$count += () = (getReverseComplementarySequence($sequence) =~ /(?=$motif)/g) if($stranded eq '' || $stranded eq 'r' || $stranded eq 'reverse');
}
print $count, "\n";

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
