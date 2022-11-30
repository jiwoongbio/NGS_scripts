#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw/sum/;

my (@fastqFileList) = @ARGV;
my %baseCountHash = ();
my %positionBaseCountHash = ();
foreach my $fastqFile (@fastqFileList) {
	open(my $reader, ($fastqFile =~ /\.gz$/ ? "gzip -dc $fastqFile |" : $fastqFile));
	while(my @lineList = getLineList($reader, 4)) {
		my $length = scalar(my @baseList = split(//, $lineList[1]));
		foreach my $index (0 .. $#baseList) {
			my $base = $baseList[$index];
			$baseCountHash{$base} += 1;
			$positionBaseCountHash{$index + 1}->{$base} += 1;
			$positionBaseCountHash{$index - $length}->{$base} += 1;
		}
	}
	close($reader);
}
my @baseList = sort keys %baseCountHash;
print join("\t", 'position', @baseList), "\n";
foreach my $position (sort {$a <=> $b} grep {$_ > 0} keys %positionBaseCountHash) {
	my %baseCountHash = %{$positionBaseCountHash{$position}};
	my $sum = sum(my @countList = map {defined($_) ? $_ : 0} @baseCountHash{@baseList});
	print join("\t", $position, map {$_ / $sum} @countList), "\n";
}
foreach my $position (sort {$a <=> $b} grep {$_ < 0} keys %positionBaseCountHash) {
	my %baseCountHash = %{$positionBaseCountHash{$position}};
	my $sum = sum(my @countList = map {defined($_) ? $_ : 0} @baseCountHash{@baseList});
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
