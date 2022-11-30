#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($inputFile) = @ARGV;
open(my $reader, $inputFile);
while(my $line = <$reader>) {
	chomp($line);
	my ($taxonNames, $taxonIds, $count) = split(/\t/, $line, -1);
	my @taxonNameList = split(/;/, $taxonNames);
	my @taxonIdList = split(/;/, $taxonIds);
	foreach my $index (0 .. $#taxonIdList) {
		if($taxonIdList[$index] ne '') {
			print join("\t", join('', map {"$_;"} @taxonNameList[0 .. $index]), join('', map {"$_;"} @taxonIdList[0 .. $index]), $count), "\n";
		}
	}
}
close($reader);
