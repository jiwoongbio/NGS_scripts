#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'c=i' => \(my $countIndex = ''),
);
my (@regionFileList) = @ARGV;
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2n");
foreach my $regionFile (@regionFileList) {
	open(my $reader, ($regionFile =~ /\.gz$/) ? "gzip -dc $regionFile |" : $regionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my ($chromosome, $start, $end) = @tokenList;
		my $count = $countIndex eq '' ? 1 : $tokenList[$countIndex];
		print $writer join("\t", $chromosome, $start, $count), "\n";
		print $writer join("\t", $chromosome, $end + 1, -$count), "\n";
	}
	close($reader);
}
close($writer);
{
	my ($chromosome, $start, $startCount, $position, $positionCount) = ('', 0, 0, 0, 0);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		if($tokenList[0] ne $chromosome) {
			die if($positionCount != 0);
			if($positionCount != $startCount) {
				print join("\t", $chromosome, $start, $position - 1, $startCount), "\n" if($startCount != 0);
				$startCount = $positionCount;
			}
			($chromosome, $start, $position) = ($tokenList[0], 0, $tokenList[1]);
		}
		if($tokenList[1] != $position) {
			die if($positionCount < 0);
			if($positionCount != $startCount) {
				print join("\t", $chromosome, $start, $position - 1, $startCount), "\n" if($startCount != 0);
				($start, $startCount) = ($position, $positionCount);
			}
			$position = $tokenList[1];
		}
		$positionCount += $tokenList[2];
	}
	die if($positionCount != 0);
	if($positionCount != $startCount) {
		print join("\t", $chromosome, $start, $position - 1, $startCount), "\n" if($startCount != 0);
		$startCount = $positionCount;
	}
	($chromosome, $start, $position) = ('', 0, 0);
}
close($reader);
waitpid($pid, 0);
