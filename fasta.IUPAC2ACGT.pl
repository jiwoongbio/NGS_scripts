#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my %iupac2acgtHash = (
	'A' => 'A',
	'C' => 'C',
	'G' => 'G',
	'T' => 'T',
	'U' => 'T',
	'R' => 'AG',
	'Y' => 'CT',
	'S' => 'GC',
	'W' => 'AT',
	'K' => 'GT',
	'M' => 'AC',
	'B' => 'CGT',
	'D' => 'AGT',
	'H' => 'ACT',
	'V' => 'ACG',
	'N' => 'ACGT',
	'.' => '-',
	'-' => '-',
);

my ($fastaFile) = @ARGV;
my @nameSequenceList = ();
{
	open(my $reader, ($fastaFile =~ /\.gz$/ ? "gzip -dc $fastaFile |" : $fastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^>//) {
			push(@nameSequenceList, [$line, '']);
		} else {
			$nameSequenceList[-1]->[1] .= $line;
		}
	}
	close($reader);
}
foreach(@nameSequenceList) {
	my ($name, $sequence) = @$_;
	my @sequenceList = ('');
	foreach my $index (0 .. length($sequence) - 1) {
		my @nextSequenceList = ();
		foreach my $base (split(//, $iupac2acgtHash{substr($sequence, $index, 1)})) {
			push(@nextSequenceList, map {"$_$base"} @sequenceList);
		}
		@sequenceList = @nextSequenceList;
	}
	foreach my $index (0 .. $#sequenceList) {
		if(scalar(@sequenceList) > 1) {
			printf(">%s_%s\n", $name, $index + 1);
		} else {
			print ">$name\n";
		}
		print "$sequenceList[$index]\n";
	}
}
