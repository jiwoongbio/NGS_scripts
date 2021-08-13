#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min);

my ($phylipOutputFile, @nameList) = @ARGV;
my %nameHash = ();
foreach my $number (1 .. scalar(@nameList)) {
	my $name = $nameList[$number - 1];
	$nameHash{sprintf('s%08d', $number)} = $name;
}
open(my $reader, $phylipOutputFile);
while(my $line = <$reader>) {
	chomp($line);
	$line =~ s/(s[0-9]{8})/$nameHash{$1}/g;
	print $line, "\n";
}
close($reader);
