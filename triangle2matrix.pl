#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($triangleFile) = @ARGV;
my @nameList = ();
my @distanceListList = ();
my $index = 0;
open(my $reader, $triangleFile);
while(my $line = <$reader>) {
	chomp($line);
	my ($name, @distanceList) = split(/\t/, $line);
	push(@nameList, $name);
	foreach(0 .. $#distanceList) {
		next if($distanceList[$_] eq '');
		$distanceListList[$index]->[$_] += $distanceList[$_];
		$distanceListList[$_]->[$index] += $distanceList[$_];
	}
	$index += 1;
}
close($reader);

print join("\t", '', @nameList), "\n";
foreach my $index (0 .. $#nameList) {
	print join("\t", $nameList[$index], map {defined($_) ? $_ : 0} @{$distanceListList[$index]}[0 .. $#nameList]), "\n";
}
