#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use List::Util qw(all min max);

my (@regionFileList) = @ARGV;
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2n -k3,3n");
foreach my $regionFile (grep {-s $_} grep {-r $_} @regionFileList) {
	open(my $reader, "sort -t '\t' -k1,1 -k2,2n -k3,3n $regionFile |");
	chomp(my $line = <$reader>);
	my ($chromosome, $start, $end) = split(/\t/, $line, -1);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		if($tokenList[0] ne $chromosome || $end < $tokenList[1]) {
			print $writer join("\t", $chromosome, $start, $end, $regionFile), "\n";
			($chromosome, $start, $end) = @tokenList;
		} elsif($end < $tokenList[2]) {
			$end = $tokenList[2];
		}
	}
	print $writer join("\t", $chromosome, $start, $end, $regionFile), "\n";
	close($reader);
}
close($writer);
my %regionFileChromosomeStartEndHash = ();
while(my $line = <$reader>) {
	chomp($line);
	my ($chromosome, $start, $end, $regionFile) = split(/\t/, $line, -1);
	$regionFileChromosomeStartEndHash{$regionFile} = [$chromosome, $start, $end];
	if(all {defined($_) && $_->[0] eq $chromosome} (my @chromosomeStartEndList = @regionFileChromosomeStartEndHash{@regionFileList})) {
		my ($start, $end) = (max(map {$_->[1]} @chromosomeStartEndList), min(map {$_->[2]} @chromosomeStartEndList));
		if($start <= $end) {
			print join("\t", $chromosome, $start, $end), "\n";
		}
	}
}
close($reader);
waitpid($pid, 0);
