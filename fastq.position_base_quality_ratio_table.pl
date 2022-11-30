#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw/sum/;

my (@fastqFileList) = @ARGV;
my %baseQualityCountHash = ();
my %positionBaseQualityCountHash = ();
foreach my $fastqFile (@fastqFileList) {
	open(my $reader, ($fastqFile =~ /\.gz$/ ? "gzip -dc $fastqFile |" : $fastqFile));
	while(my @lineList = getLineList($reader, 4)) {
		my $length = scalar(my @baseQualityList = split(//, $lineList[3]));
		foreach my $index (0 .. $#baseQualityList) {
			my $baseQuality = ord($baseQualityList[$index]);
			$baseQualityCountHash{$baseQuality} += 1;
			$positionBaseQualityCountHash{$index + 1}->{$baseQuality} += 1;
			$positionBaseQualityCountHash{$index - $length}->{$baseQuality} += 1;
		}
	}
	close($reader);
}
my @baseQualityList = sort {$a <=> $b} keys %baseQualityCountHash;
print join("\t", 'position', map {"Q$_"} @baseQualityList), "\n";
foreach my $position (sort {$a <=> $b} grep {$_ > 0} keys %positionBaseQualityCountHash) {
	my %baseQualityCountHash = %{$positionBaseQualityCountHash{$position}};
	my $sum = sum(my @countList = map {defined($_) ? $_ : 0} @baseQualityCountHash{@baseQualityList});
	print join("\t", $position, map {$_ / $sum} @countList), "\n";
}
foreach my $position (sort {$a <=> $b} grep {$_ < 0} keys %positionBaseQualityCountHash) {
	my %baseQualityCountHash = %{$positionBaseQualityCountHash{$position}};
	my $sum = sum(my @countList = map {defined($_) ? $_ : 0} @baseQualityCountHash{@baseQualityList});
	print join("\t", $position, map {$_ / $sum} @countList), "\n";
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
