#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'v' => \(my $invertMatch = ''),
	'h' => \(my $header = ''),
	'i' => \(my $ignoreCase = ''),
);
my ($tableFile1, $indexes1, $tableFile2, $indexes2) = @ARGV;
my @indexList1 = eval($indexes1);
my @indexList2 = eval($indexes2);
my %searchHash = ();
{
	open(my $reader, $tableFile1);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = $line eq '' ? ('') : split(/\t/, $line, -1);
		push(@tokenList, '') if(scalar(@tokenList) == 0);
		my $search = join("\t", @tokenList[@indexList1]);
		$searchHash{$ignoreCase ? uc($search) : $search} = 1;
	}
	close($reader);
}
{
	open(my $reader, $tableFile2);
	if($header) {
		chomp(my $line = <$reader>);
		print "$line\n";
	}
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = $line eq '' ? ('') : split(/\t/, $line, -1);
		push(@tokenList, '') if(scalar(@tokenList) == 0);
		my $search = join("\t", @tokenList[@indexList2]);
		if(defined($searchHash{$ignoreCase ? uc($search) : $search})) {
			print "$line\n" unless($invertMatch);
		} else {
			print "$line\n" if($invertMatch);
		}
	}
	close($reader);
}
