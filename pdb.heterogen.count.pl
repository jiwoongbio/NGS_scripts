#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);

my ($pdbFile) = @ARGV;
my %residueHash = ();
open(my $reader, ($pdbFile =~ /\.gz$/ ? "gzip -dc $pdbFile |" : $pdbFile));
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^HETATM/) {
		my $residue = substr($line, 17, 10);
		$residueHash{$residue} = 1;
	}
}
close($reader);

my %residueNameCountHash = ();
foreach my $residue (sort keys %residueHash) {
	my $residueName = substr($residue, 0, 3);
	$residueNameCountHash{$residueName} += 1;
}

foreach my $residueName (sort keys %residueNameCountHash) {
	print join("\t", $residueName, $residueNameCountHash{$residueName}), "\n";
}
