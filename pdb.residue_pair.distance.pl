#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);

my @aaResidueNameList = ('ALA', 'ASX', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR', 'GLX');
my %aaResidueNameHash = map {$_ => 1} @aaResidueNameList;

my @proteinBackboneAtomNameList = (' N  ', ' CA ', ' C  ', ' O  ');
my %proteinBackboneAtomNameHash = map {$_ => 1} @proteinBackboneAtomNameList;

my ($pdbFile) = @ARGV;
my $model = '';
my %modelResidueXYZListHash = ();
open(my $reader, ($pdbFile =~ /\.gz$/ ? "gzip -dc $pdbFile |" : $pdbFile));
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^MODEL /) {
		$model = substr($line, 10, 4);
	}
	if($line =~ /^ATOM  /) {
		my $elementSymbol = substr($line, 76, 2);
		if($elementSymbol ne ' H') {
			my $residueName = substr($line, 17, 3);
			my $atomName = substr($line, 12, 4);
			if(!($aaResidueNameHash{$residueName} && $proteinBackboneAtomNameHash{$atomName}) || ($residueName eq 'GLY' && $atomName eq ' CA ')) {
				my $residue = substr($line, 17, 10);
				my $x = substr($line, 30, 8) + 0;
				my $y = substr($line, 38, 8) + 0;
				my $z = substr($line, 46, 8) + 0;
				push(@{$modelResidueXYZListHash{$model}->{$residue}}, [$x, $y, $z]);
			}
		}
	}
	if($line =~ /^HETATM/) {
		my $elementSymbol = substr($line, 76, 2);
		if($elementSymbol ne ' H') {
			my $residueName = substr($line, 17, 3);
			if($residueName ne 'HOH') {
				my $residue = substr($line, 17, 10);
				my $x = substr($line, 30, 8) + 0;
				my $y = substr($line, 38, 8) + 0;
				my $z = substr($line, 46, 8) + 0;
				push(@{$modelResidueXYZListHash{$model}->{$residue}}, [$x, $y, $z]);
			}
		}
	}
}
close($reader);

my %residuePairDistanceHash = ();
foreach my $model (sort keys %modelResidueXYZListHash) {
	my %residueXYZListHash = %{$modelResidueXYZListHash{$model}};
	my @residueList = sort keys %residueXYZListHash;
	foreach my $residue1 (@residueList) {
		foreach my $residue2 (@residueList) {
			if($residue1 ne $residue2) {
				my $residuePair = join("\t", $residue1, $residue2);
				foreach my $xyz1 (@{$residueXYZListHash{$residue1}}) {
					foreach my $xyz2 (@{$residueXYZListHash{$residue2}}) {
						my $distance = sum(0, map {($xyz1->[$_] - $xyz2->[$_]) ** 2} 0 .. 2) ** 0.5;
						$residuePairDistanceHash{$residuePair} = $distance if(!defined($residuePairDistanceHash{$residuePair}) || $distance < $residuePairDistanceHash{$residuePair});
					}
				}
			}
		}
	}
}

foreach my $residuePair (sort keys %residuePairDistanceHash) {
	print join("\t", $residuePair, $residuePairDistanceHash{$residuePair}), "\n";
}
