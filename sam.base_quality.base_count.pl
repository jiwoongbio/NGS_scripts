#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my (@samFileList) = @ARGV;
my %baseQualityCountHash = ();
foreach my $samFile (@samFileList) {
	open(my $reader, ($samFile =~ /\.gz$/ ? "gzip -dc $samFile |" : $samFile));
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^@/);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		$baseQualityCountHash{$_} += 1 foreach(split(//, $tokenHash{'qual'}));
	}
	close($reader);
}
%baseQualityCountHash = map {ord($_) => $baseQualityCountHash{$_}} keys %baseQualityCountHash;
my @baseQualityList = sort {$a <=> $b} keys %baseQualityCountHash;
foreach my $baseQuality ($baseQualityList[0] .. $baseQualityList[-1]) {
	my $count = $baseQualityCountHash{$baseQuality};
	print join("\t", chr($baseQuality), $baseQuality, defined($count) ? $count : 0), "\n";
}
