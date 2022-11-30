#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(max);

my ($tableFile) = @ARGV;
my @tokenListList = ();
open(my $reader, $tableFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	push(@tokenList, '') if(scalar(@tokenList) == 0);
	push(@tokenListList, \@tokenList);
}
close($reader);
foreach my $index (0 .. max(map {$#$_} @tokenListList)) {
	print join("\t", map {defined($_) ? $_ : ''} map {$_->[$index]} @tokenListList), "\n";
}
