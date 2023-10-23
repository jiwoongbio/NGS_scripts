#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;

my ($mutationRateMatrixFile) = @ARGV;
open(my $reader, $mutationRateMatrixFile);
chomp(my $line = <$reader>);
my @tokenList = split(/\t/, $line, -1);
print join("\t", @tokenList), "\n";
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	print join("\t", $tokenList[0], map {jukesCantorDistance($_)} @tokenList[1 .. $#tokenList]), "\n";
}
close($reader);

sub jukesCantorDistance {
	my ($mutationRate) = @_;
	return (-3 / 4) * log(1 - (4 / 3) * $mutationRate);
}
