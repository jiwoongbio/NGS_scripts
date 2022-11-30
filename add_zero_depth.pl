#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($depthFile) = @ARGV;
my ($previousChromosome, $previousPosition);
open(my $reader, ($depthFile =~ /\.gz$/ ? "gzip -dc $depthFile |" : $depthFile));
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome, $position, $depth) = split(/\t/, $line, -1);
	if(defined($previousChromosome) && $previousChromosome eq $chromosome) {
		my $depth = 0;
		print join("\t", $chromosome, $_, 0), "\n" foreach($previousPosition + 1 .. $position - 1);
	}
	print join("\t", $chromosome, $position, $depth), "\n";
	($previousChromosome, $previousPosition) = ($chromosome, $position);
}
close($reader);
