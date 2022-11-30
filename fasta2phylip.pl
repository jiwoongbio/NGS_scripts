#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min);

my ($fastaFile) = @ARGV;
my @nameSequenceList = ();
{
	open(my $reader, ($fastaFile =~ /\.gz$/ ? "gzip -dc $fastaFile |" : $fastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			push(@nameSequenceList, [$1, '']);
		} else {
			$nameSequenceList[-1]->[1] .= $line;
		}
	}
	close($reader);
}
my $number = scalar(@nameSequenceList);
my $length = min(map {length($_->[1])} @nameSequenceList);
print join(' ', $number, $length), "\n";
foreach my $number (1 .. $number) {
	my ($name, $sequence) = @{$nameSequenceList[$number - 1]};
	print join(' ', sprintf('s%08d', $number), $sequence), "\n";
}
