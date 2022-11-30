#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
);
my (@samFileList) = @ARGV;
my %chromosomeMatchDepthHash = ();
foreach my $samFile (@samFileList) {
	open(my $reader, "samtools view -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $samFile |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		$chromosomeMatchDepthHash{$tokenHash{'rname'}} += sum(0, $tokenHash{'MD:Z'} =~ /([0-9]+)/g) if(defined($tokenHash{'MD:Z'}));
	}
	close($reader);
}
print join("\t", $_, $chromosomeMatchDepthHash{$_}), "\n" foreach(sort keys %chromosomeMatchDepthHash);
