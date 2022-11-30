#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);

my ($tableFile) = @ARGV;
open(my $reader, ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile |" : $tableFile));
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line, -1);
my @sampleList = @columnList[1 .. $#columnList];
my @valueListList = ();
my @indexList = ();
my $index = 0;
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my @valueList = @tokenList[1 .. $#tokenList];
	@valueList = map {$_ eq '' ? 0 : $_} @valueList;
	push(@{$valueListList[$_]}, $valueList[$_]) foreach(0 .. $#valueList);
	push(@indexList, $index);
	$index += 1;
}
close($reader);
print join("\t", '', @sampleList), "\n";
foreach my $index1 (0 .. $#sampleList) {
	my $valueList1 = $valueListList[$index1];
	my @distanceList = ();
	foreach my $index2 (0 .. $#sampleList) {
		my $valueList2 = $valueListList[$index2];
		my $distance = sum(0, map {($valueList1->[$_] - $valueList2->[$_]) ** 2} @indexList) ** 0.5;
		push(@distanceList, $distance);
	}
	print join("\t", $sampleList[$index1], @distanceList), "\n";
}
