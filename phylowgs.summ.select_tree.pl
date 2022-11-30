#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;

my ($summFile, @treeIndexList) = @ARGV;
my $summ = '';
{
	open(my $reader, $summFile);
	$summ .= $_ while(<$reader>);
	close($reader);
}
$summ = decode_json($summ);
my $trees;
foreach my $index (@treeIndexList) {
	$trees->{$index} = $summ->{'trees'}->{$index};
}
$summ->{'trees'} = $trees;
print encode_json($summ);
