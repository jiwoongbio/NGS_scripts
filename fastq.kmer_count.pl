#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'q=i' => \(my $minimumBaseQuality = 0),
);
my ($k, @fastqFileList) = @ARGV;
my %kmerCountHash = ();
my $count = 0;
foreach my $fastqFile (@fastqFileList) {
	open(my $reader, ($fastqFile =~ /\.gz$/ ? "gzip -dc $fastqFile |" : $fastqFile));
	while(my @lineList = getLineList($reader, 4)) {
		my $sequence = $lineList[1];
		my @baseQualityList = map {ord($_)} split(//, $lineList[3]);
		my %indexHash = ();
		foreach my $index (0 .. $#baseQualityList) {
			if($baseQualityList[$index] < $minimumBaseQuality) {
				$indexHash{$_} = 1 foreach($index - $k + 1 .. $index + $k - 1);
			}
		}
		foreach my $index (length($sequence) - $k) {
			next if($indexHash{$index});
			my $kmer = substr($sequence, $index, $k);
			$kmerCountHash{$kmer} += 1;
			$count += 1;
		}
	}
	close($reader);
}
foreach my $kmer (sort keys %kmerCountHash) {
	print join("\t", $kmer, $kmerCountHash{$kmer}, $kmerCountHash{$kmer} / $count), "\n";
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
