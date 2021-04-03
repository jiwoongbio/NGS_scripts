#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	's' => \(my $sorted = ''),
	'g=i' => \(my $gap = 0),
);
my ($regionFile) = @ARGV;
open(my $reader, $sorted ? $regionFile : "sort --field-separator=\$'\\t' -k4 -k1,1 -k2,2n -k3,3n $regionFile |");
my @previousTokenListList = ();
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, 4);
	if(@previousTokenListList) {
		my @remainedPreviousTokenListList = ();
		foreach my $previousTokenList (@previousTokenListList) {
			if($previousTokenList->[0] ne $tokenList[0] || $previousTokenList->[2] + $gap + 1 < $tokenList[1]) {
				print join("\t", @$previousTokenList), "\n";
			} else {
				push(@remainedPreviousTokenListList, $previousTokenList);
			}
		}
		@previousTokenListList = @remainedPreviousTokenListList;
	}
	my $added = '';
	foreach my $previousTokenList (@previousTokenListList) {
		if((defined($previousTokenList->[3]) && defined($tokenList[3]) && $previousTokenList->[3] eq $tokenList[3]) || (!defined($previousTokenList->[3]) && !defined($tokenList[3]))) {
			$previousTokenList->[2] = $tokenList[2] if($previousTokenList->[2] < $tokenList[2]);
			$added = 1;
		}
	}
	push(@previousTokenListList, \@tokenList) if($added eq '');
}
print join("\t", @$_), "\n" foreach(@previousTokenListList);
close($reader);
