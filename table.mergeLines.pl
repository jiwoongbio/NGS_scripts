#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum max min);
use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
	's' => \(my $sorted = ''),
	'u' => \(my $unique = ''),
	'm' => \(my $multiple = ''),
	'f=s' => \(my $function = 'join'),
	'd=s' => \(my $delimiter = ','),
	'N=i' => \(my $number = ''),
	'D=s' => \(my $default = ''),
);
$function = 'lower_quartile' if($function eq 'Q1');
$function = 'median' if($function eq 'Q2');
$function = 'upper_quartile' if($function eq 'Q3');

if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] table.txt indexes [key_index value_index [key1 ...]] > table.line_merged.txt

Options: -h       display this help message
         -s       input table is sorted
         -u       output unique values
         -m       output multiple values
         -f STR   merge function [$function]
         -d STR   delimiter [$delimiter]
         -N INT   number of values
         -D STR   default value

EOF
}
my ($tableFile, $indexes, $keyIndex, $valueIndex, @keyList) = @ARGV;
my %indexHash = map {$_ => 1} (my @indexList = eval($indexes));
if(defined($valueIndex) && scalar(@keyList) == 0) {
	my %keyHash = ();
	open(my $reader, ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile |" : $tableFile));
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
	open(my $reader, $sorted ? ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile |" : $tableFile) : ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile | sort -t '\t' @sortKeyOptionList |" : "sort -t '\t' @sortKeyOptionList $tableFile |"));
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
	if($number eq '') {
	} elsif(scalar(@tokenList) > $number) {
		@tokenList = @tokenList[0 .. $number - 1];
	} elsif(scalar(@tokenList) < $number) {
		push(@tokenList, ($default) x ($number - scalar(@tokenList)));
	}
	my @uniqueTokenList = unique(@tokenList);
	@tokenList = @uniqueTokenList if($unique);
	if($function eq 'join') {
		if($multiple eq '' && scalar(@uniqueTokenList) == 1) {
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
	@sortedTokenList = @sortedTokenList[0, grep {$sortedTokenList[$_ - 1] ne $sortedTokenList[$_]} 1 .. $#sortedTokenList] if(@sortedTokenList);
	return @sortedTokenList;
}

sub mean {
	my @tokenList = @_;
	return sum(@tokenList) / scalar(@tokenList);
}

sub sd {
	my @tokenList = @_;
	my $mean = mean(@tokenList);
	return (sum(map {($_ - $mean) ** 2} @tokenList) / scalar(@tokenList)) ** 0.5;
}

sub medianIndexList {
	my ($length) = @_;
	if($length % 2 == 0) {
		return ($length / 2 - 1, $length / 2);
	} else {
		return (($length - 1) / 2);
	}
}

sub median {
	my @tokenList = @_;
	@tokenList = sort {$a <=> $b} @tokenList;
	return mean(@tokenList[medianIndexList(scalar(@tokenList))]);
}

sub lower_quartile {
	my @tokenList = @_;
	@tokenList = sort {$a <=> $b} @tokenList;
	my @medianIndexList = medianIndexList(scalar(@tokenList));
	@tokenList = @tokenList[0 .. ($medianIndexList[-1] - 1)];
	return mean(@tokenList[medianIndexList(scalar(@tokenList))]);
}

sub upper_quartile {
	my @tokenList = @_;
	@tokenList = sort {$a <=> $b} @tokenList;
	my @medianIndexList = medianIndexList(scalar(@tokenList));
	@tokenList = @tokenList[($medianIndexList[0] + 1) .. $#tokenList];
	return mean(@tokenList[medianIndexList(scalar(@tokenList))]);
}
