#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum max min);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	's' => \(my $sorted = ''),
	'u' => \(my $unique = ''),
	'f=s' => \(my $function = 'join'),
	'd=s' => \(my $delimiter = ','),
);
my ($tableFile, $indexes, $keyIndex, $valueIndex, @keyList) = @ARGV;
my %indexHash = map {$_ => 1} (my @indexList = eval($indexes));
if(defined($valueIndex) && scalar(@keyList) == 0) {
	my %keyHash = ();
	open(my $reader, $tableFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		$keyHash{$tokenList[$keyIndex]} = 1;
	}
	close($reader);
	@keyList = sort keys %keyHash;
}
{
	my @sortKeyOptionList = map {sprintf('-k%d,%d', $_ + 1, $_ + 1)} @indexList;
	open(my $reader, $sorted ? $tableFile : "sort -t '\t' @sortKeyOptionList $tableFile |");
	my ($previousKey, @tokenListList) = ('');
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		if((my $key = join("\t", @tokenList[@indexList])) ne $previousKey) {
			printLine(@tokenListList) if(@tokenListList);
			($previousKey, @tokenListList) = ($key);
		}
		push(@tokenListList, \@tokenList);
	}
	printLine(@tokenListList) if(@tokenListList);
	close($reader);
}

sub printLine {
	my (@tokenListList) = @_;
	my $length = max(map {scalar(@$_)} @tokenListList);
	my %tokenHash = ();
	if(@keyList) {
		foreach(@tokenListList) {
			$tokenHash{$_->[$keyIndex]} = $_->[$valueIndex];
			($_->[$keyIndex], $_->[$valueIndex]) = ('', '');
		}
		@tokenListList = map {[split(/\t/, $_, -1)]} unique(map {join("\t", @$_)} @tokenListList);
	}
	my @tokenList = ();
	for(my $index = 0; $index < $length; $index++) {
		if($indexHash{$index}) {
			$tokenList[$index] = $tokenListList[0]->[$index];
		} else {
			$tokenList[$index] = getToken(grep {defined} map {$_->[$index]} @tokenListList);
		}
	}
	if(@keyList) {
		$tokenList[$valueIndex] = join("\t", map {defined($_) ? $_ : ''} @tokenHash{@keyList});
		splice(@tokenList, $keyIndex, 1);
	}
	print join("\t", @tokenList), "\n";
}

sub getToken {
	my @tokenList = @_;
	my @uniqueTokenList = unique(@tokenList);
	@tokenList = @uniqueTokenList if($unique);
	if($function eq 'join') {
		if(scalar(@uniqueTokenList) == 1) {
			return $uniqueTokenList[0];
		} else {
			return join($delimiter, @tokenList);
		}
	} else {
		return eval("$function(\@tokenList)");
	}
}

sub unique {
	my @sortedTokenList = sort @_;
	@sortedTokenList = @sortedTokenList[0, grep {$sortedTokenList[$_ - 1] ne $sortedTokenList[$_]} 1 .. $#sortedTokenList];
	return @sortedTokenList;
}

sub mean {
	my @tokenList = @_;
	return sum(@tokenList) / scalar(@tokenList);
}

sub median {
	my @tokenList = @_;
	@tokenList = sort {$a <=> $b} @tokenList;
	if((my $length = scalar(@tokenList)) % 2) {
		return ($tokenList[$length / 2 - 1] + $tokenList[$length / 2]) / 2;
	} else {
		return $tokenList[($length - 1) / 2];
	}
}
