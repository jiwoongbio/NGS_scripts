#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my %aaResidueNameOneLetterHash = ('ALA' => 'A', 'ASX' => 'B', 'CYS' => 'C', 'ASP' => 'D', 'GLU' => 'E', 'PHE' => 'F', 'GLY' => 'G', 'HIS' => 'H', 'ILE' => 'I', 'LYS' => 'K', 'LEU' => 'L', 'MET' => 'M', 'ASN' => 'N', 'PRO' => 'P', 'GLN' => 'Q', 'ARG' => 'R', 'SER' => 'S', 'THR' => 'T', 'VAL' => 'V', 'TRP' => 'W', 'TYR' => 'Y', 'GLX' => 'Z');

GetOptions(
	'fastaLineLength=i' => \(my $fastaLineLength = 80),
);
my ($pdbFile, @chainList) = @ARGV;
my %chainHash = ();
foreach my $chain (@chainList) {
	$chainHash{$chain} = 1;
}
my %chainResidueNumberNameHash = ();
open(my $reader, ($pdbFile =~ /\.gz$/ ? "gzip -dc $pdbFile |" : $pdbFile));
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^ATOM  /) {
		my $chain = substr($line, 21, 1);
		if(scalar(@chainList) == 0 || $chainHash{$chain}) {
			my $insertionResidueCode = substr($line, 26, 1);
			if($insertionResidueCode eq ' ') {
				my $residueNumber = int(substr($line, 22, 4));
				my $residueName = substr($line, 17, 3);
				$chainResidueNumberNameHash{$chain}->{$residueNumber} = $residueName unless(defined($chainResidueNumberNameHash{$chain}->{$residueNumber}));
			}
		}
	}
}
close($reader);

foreach my $chain (sort keys %chainResidueNumberNameHash) {
	my $residueNumberNameHash = $chainResidueNumberNameHash{$chain};
	my @residueNumberList = sort {$a <=> $b} keys %$residueNumberNameHash;
	my ($startResidueNumber, $endResidueNumber) = @residueNumberList[0, -1];
	$startResidueNumber = 1 if($startResidueNumber < 1);
	my $sequence = '';
	foreach my $residueNumber ($startResidueNumber .. $endResidueNumber) {
		if(defined(my $residueName = $residueNumberNameHash->{$residueNumber})) {
			if(defined(my $residueOneLetter = $aaResidueNameOneLetterHash{$residueName})) {
				$sequence .= $residueOneLetter;
			} else {
				$sequence .= 'X';
			}
		} else {
			$sequence .= 'X';
		}
	}
	$chain =~ s/ //g;
	my $offset = $startResidueNumber - 1;
	print ">$chain|$offset\n";
	for(my $index = 0; $index < length($sequence); $index += $fastaLineLength) {
		print substr($sequence, $index, $fastaLineLength), "\n";
	}
}
