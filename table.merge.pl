#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum max);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'u' => \(my $unique = ''),
	'f=s' => \(my $function = 'join'),
	'd=s' => \(my $delimiter = ','),
);
my (@regionFileList) = @ARGV;
my @readerList = ();
foreach my $regionFile (@regionFileList) {
	open(my $reader, $regionFile);
	push(@readerList, $reader);
}
while(my @tokenListList = getTokenListList()) {
	my @tokenList = ();
	foreach my $index (0 .. max(map {$#$_} @tokenListList)) {
		$tokenList[$index] = getToken(grep {defined} map {$_->[$index]} @tokenListList);
	}
	print join("\t", @tokenList), "\n";
}
close($_) foreach(@regionFileList);

sub getTokenListList {
	my @tokenListList = ();
	foreach my $reader (@readerList) {
		if(defined(my $line = <$reader>)) {
			chomp($line);
			my @tokenList = split(/\t/, $line, -1);
			push(@tokenListList, \@tokenList);
		}
	}
	return @tokenListList;
}

sub getToken {
	my @tokenList = @_;
	my @uniqueTokenList = unique(@tokenList);
	return $uniqueTokenList[0] if(scalar(@uniqueTokenList) == 1);
	@tokenList = @uniqueTokenList if($unique);
	if($function eq 'join') {
		return join($delimiter, @tokenList);
	} else {
		return eval("$function(\@tokenList)");
	}
}

sub unique {
	my @sorted = sort @_;
	return @sorted[0, grep {$sorted[$_ - 1] ne $sorted[$_]} 1 .. $#sorted];
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
