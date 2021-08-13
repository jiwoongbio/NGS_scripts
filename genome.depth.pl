#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use bignum;
use Getopt::Long qw(:config no_ignore_case);

my @coverageDepthList = ();
GetOptions(
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'c=i' => \@coverageDepthList,
);
my ($bamFile) = @ARGV;
my $totalLength = 0;
{
	open(my $reader, "samtools view -H $bamFile |");
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		if($tokenList[0] eq '@SQ') {
			my %tokenHash = map {$_->[0] => $_->[1]} map {[split(/:/, $_, 2)]} @tokenList[1 .. $#tokenList];
			$totalLength += $tokenHash{'LN'};
		}
	}
	close($reader);
}
my $totalDepth = 0;
my @totalCoverageList = (0) x scalar(@coverageDepthList);
{
	open(my $reader, "samtools view -u -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile | samtools depth - |");
	while(my $line = <$reader>) {
		chomp($line);
		my ($chromosome, $position, $positionDepth) = split(/\t/, $line, -1);
		$totalDepth += $positionDepth;
		$totalCoverageList[$_] += 1 foreach(grep {$positionDepth >= $coverageDepthList[$_]} 0 .. $#coverageDepthList);
	}
	close($reader);
}
print join("\t", 'Mean depth', $totalDepth / $totalLength), "\n";
print join("\t", "% Depth >= $coverageDepthList[$_]", sprintf('%.1f%%', $totalCoverageList[$_] / $totalLength * 100)), "\n" foreach(0 .. $#totalCoverageList);
