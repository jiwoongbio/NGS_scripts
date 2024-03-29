#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
	'i' => \(my $ignoreCase = ''),
	'm' => \(my $matchedOnly = ''),
	'o' => \(my $once = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] table.txt key_indexes table1.txt key_indexes1 value_indexes1 [...] > table.column_added.txt

Options: -h       display this help message
         -i       ignore case
         -m       matched only
         -o       once

EOF
}
my ($tableFile, $keyIndexes, @additionTableFileKeyValueIndexesList) = @ARGV;
die "$tableFile is not readable." unless($tableFile eq '-' || $tableFile =~ /\|\s*$/ || -r $tableFile);
die unless(scalar(@additionTableFileKeyValueIndexesList) % 3 == 0);
my (@additionListHashList, @spacerList) = ();
for(my $index = 0; $index < scalar(@additionTableFileKeyValueIndexesList); $index += 3) {
	my ($additionTableFile, $additionKeyIndexes, $additionValueIndexes) = @additionTableFileKeyValueIndexesList[$index .. $index + 2];
	die "$additionTableFile is not readable." unless($additionTableFile eq '-' || $additionTableFile =~ /\|\s*$/ || -r $additionTableFile);
	my @additionKeyIndexList = eval($additionKeyIndexes);
	my @additionValueIndexList = eval($additionValueIndexes);
	my %additionListHash = ();
	open(my $reader, $additionTableFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = $line eq '' ? ('') : split(/\t/, $line, -1);
		my $key = join("\t", @tokenList[@additionKeyIndexList]);
		push(@{$additionListHash{$ignoreCase ? uc($key) : $key}}, join("\t", @tokenList[@additionValueIndexList]));
	}
	close($reader);
	push(@additionListHashList, \%additionListHash);
	push(@spacerList, join("\t", ('') x scalar(@additionValueIndexList)));
}
my @keyIndexList = eval($keyIndexes);
open(my $reader, $tableFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = $line eq '' ? ('') : split(/\t/, $line, -1);
	my $key = join("\t", @tokenList[@keyIndexList]);
	my @additionListList = map {$_->{$ignoreCase ? uc($key) : $key}} @additionListHashList;
	if($once) {
		delete $_->{$ignoreCase ? uc($key) : $key} foreach(@additionListHashList);
	}
	next if(scalar(my @unmatchedIndexList = grep {!defined($additionListList[$_])} 0 .. $#additionListList) > 0 && $matchedOnly);
	@additionListList[@unmatchedIndexList] = map {[$_]} @spacerList[@unmatchedIndexList];
	my $tokenListList = [\@tokenList];
	foreach my $additionList (@additionListList) {
		my @tokenListList = ();
		foreach my $tokenList (@$tokenListList) {
			push(@tokenListList, map {[@$tokenList, $_]} @$additionList);
		}
		$tokenListList = \@tokenListList;
	}
	foreach my $tokenList (@$tokenListList) {
		print join("\t", @$tokenList), "\n";
	}
}
close($reader);
