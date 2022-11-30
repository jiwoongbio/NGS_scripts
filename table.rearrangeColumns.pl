#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my @excludeIndexesList = ();
GetOptions(
	'e=s' => \@excludeIndexesList,
	'c' => \(my $useColumnInsteadOfIndexes = ''),
	'f=s' => \(my $columnFile = ''),
);
my ($tableFile, @indexesList) = @ARGV;
if($columnFile ne '') {
	open(my $reader, $columnFile);
	while(my $line = <$reader>) {
		chomp($line);
		push(@indexesList, $line);
	}
	close($reader);
}
if(@indexesList) {
	open(my $reader, $tableFile);
	my @indexList = ();
	my %excludeIndexHash = ();
	if($useColumnInsteadOfIndexes) {
		chomp(my $line = <$reader>);
		my @tokenList = split(/\t/, $line, -1);
		my %columnIndexHash = map {$tokenList[$_] => $_} 0 .. $#tokenList;
		if((my $column = $tokenList[0]) =~ s/^#//) {
			$columnIndexHash{$column} = 0;
		}
		push(@indexList, @columnIndexHash{@indexesList});
		$excludeIndexHash{$_} = 1 foreach(@columnIndexHash{@excludeIndexesList});
		print join("\t", @tokenList[grep {!$excludeIndexHash{$_}} @indexList]), "\n";
	} else {
		push(@indexList, map {eval($_)} @indexesList);
		$excludeIndexHash{$_} = 1 foreach(map {eval($_)} @excludeIndexesList);
	}
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		print join("\t", @tokenList[grep {!$excludeIndexHash{$_}} @indexList]), "\n";
	}
	close($reader);
} else {
	open(my $reader, $tableFile);
	my %excludeIndexHash = ();
	if($useColumnInsteadOfIndexes) {
		chomp(my $line = <$reader>);
		my @tokenList = split(/\t/, $line, -1);
		my %columnIndexHash = map {$tokenList[$_] => $_} 0 .. $#tokenList;
		if((my $column = $tokenList[0]) =~ s/^#//) {
			$columnIndexHash{$column} = 0;
		}
		$excludeIndexHash{$_} = 1 foreach(@columnIndexHash{@excludeIndexesList});
		my @indexList = 0 .. $#tokenList;
		print join("\t", @tokenList[grep {!$excludeIndexHash{$_}} @indexList]), "\n";
	} else {
		$excludeIndexHash{$_} = 1 foreach(map {eval($_)} @excludeIndexesList);
	}
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my @indexList = 0 .. $#tokenList;
		print join("\t", @tokenList[grep {!$excludeIndexHash{$_}} @indexList]), "\n";
	}
	close($reader);
}
