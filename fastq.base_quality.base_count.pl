#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my (@fastqFileList) = @ARGV;
my %baseQualityCountHash = ();
foreach my $fastqFile (@fastqFileList) {
	open(my $reader, ($fastqFile =~ /\.gz$/ ? "gzip -dc $fastqFile |" : $fastqFile));
	while(my @lineList = getLineList($reader, 4)) {
		$baseQualityCountHash{$_} += 1 foreach(split(//, $lineList[3]));
	}
	close($reader);
}
%baseQualityCountHash = map {ord($_) => $baseQualityCountHash{$_}} keys %baseQualityCountHash;
my @baseQualityList = sort {$a <=> $b} keys %baseQualityCountHash;
foreach my $baseQuality ($baseQualityList[0] .. $baseQualityList[-1]) {
	my $count = $baseQualityCountHash{$baseQuality};
	print join("\t", chr($baseQuality), $baseQuality, defined($count) ? $count : 0), "\n";
}

sub getLineList {
	my ($reader, $number) = @_;
	my @lineList = ();
	foreach(1 .. $number) {
		if(defined(my $line = <$reader>)) {
			chomp($line);
			push(@lineList, $line);
		}
	}
	return @lineList;
}
