#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'i=s' => \(my $indexes = ''),
	'o' => \(my $substitutedOnly = ''),
	'f=s' => \(my $substituteFile = ''),
	'p' => \(my $partial = ''),
);
my @indexList = eval($indexes);
my ($tableFile, %substituteHash) = @ARGV;
if($substituteFile ne '') {
	open(my $reader, $substituteFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		$substituteHash{$tokenList[0]} = $tokenList[1];
	}
	close($reader);
}
open(my $reader, $tableFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my $substituted = '';
	if(@indexList) {
		foreach my $index (@indexList) {
			if($partial) {
				foreach(keys %substituteHash) {
					$substituted = 1 if($tokenList[$index] =~ s/$_/$substituteHash{$_}/);
				}
			} else {
				if(defined($substituteHash{$tokenList[$index]})) {
					$tokenList[$index] = $substituteHash{$tokenList[$index]};
					$substituted = 1;
				}
			}
		}
	} else {
		foreach my $index (0 .. $#tokenList) {
			if($partial) {
				foreach(keys %substituteHash) {
					$substituted = 1 if($tokenList[$index] =~ s/$_/$substituteHash{$_}/);
				}
			} else {
				if(defined($substituteHash{$tokenList[$index]})) {
					$tokenList[$index] = $substituteHash{$tokenList[$index]};
					$substituted = 1;
				}
			}
		}
	}
	print join("\t", @tokenList), "\n" if($substitutedOnly eq '' || $substituted);
}
close($reader);
