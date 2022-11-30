#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;

my ($summFile, $treeIndex) = @ARGV;
my $summ = '';
{
	open(my $reader, $summFile);
	$summ .= $_ while(<$reader>);
	close($reader);
}
$summ = decode_json($summ);
foreach my $vertex1 (sort {$a <=> $b} keys %{$summ->{'trees'}->{$treeIndex}->{'structure'}}) {
	foreach my $vertex2 (sort {$a <=> $b} @{$summ->{'trees'}->{$treeIndex}->{'structure'}->{$vertex1}}) {
		print join("\t", $vertex1, $vertex2), "\n";
	}
}
