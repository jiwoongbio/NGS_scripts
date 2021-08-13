#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($tableFile, $column, $sign, $value) = @ARGV;
open(my $reader, ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile |" : $tableFile));
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line, -1);
print join("\t", @columnList), "\n";
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my %tokenHash = ();
	@tokenHash{@columnList} = @tokenList;
	if($sign eq '=~' || $sign eq '!~') {
		if(eval(sprintf('$tokenHash{$column} %s /$value/', $sign))) {
			print join("\t", @tokenList), "\n";
		}
	} else {
		if(eval(sprintf('$tokenHash{$column} %s $value', $sign))) {
			print join("\t", @tokenList), "\n";
		}
	}
}
close($reader);
