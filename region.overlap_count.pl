#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use List::Util qw(all min max);

my (@regionFileList) = @ARGV;
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2n -k3,3n | uniq -c");
foreach my $regionFile (@regionFileList) {
	open(my $reader, ($regionFile =~ /\.gz$/) ? "gzip -dc $regionFile |" : $regionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($chromosome, $start, $end) = split(/\t/, $line, -1);
		print $writer join("\t", $chromosome, $start, 1), "\n";
		print $writer join("\t", $chromosome, $end + 1, -1), "\n";
	}
	close($reader);
}
close($writer);
{
	my ($chromosome, $start, $count) = ('', 0, 0);
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^ *([0-9]+) //) {
			my @tokenList = split(/\t/, $line, -1);
			if($tokenList[0] eq $chromosome) {
				print join("\t", $chromosome, $start, $tokenList[1] - 1, $count), "\n" if($start <= $tokenList[1] - 1 && $count > 0);
				$start = $tokenList[1];
				$count += $1 * $tokenList[2];
				die if($count < 0);
			} else {
				die if($count != 0);
				die if($tokenList[2] != 1);
				($chromosome, $start, $count) = (@tokenList[0, 1], $1);
			}
		}
	}
	die if($count != 0);
}
close($reader);
waitpid($pid, 0);
