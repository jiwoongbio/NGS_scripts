#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   fastq.read_length.read_count.pl [options] input.fastq [...]

Options: -h       display this help message

EOF
}
my (@fastqFileList) = @ARGV;
foreach my $fastqFile (@fastqFileList) {
	my %readLengthCountHash = ();
	my $error = '';
	open(my $reader, ($fastqFile =~ /\.gz$/ ? "gzip -dc $fastqFile |" : $fastqFile));
	while(my @lineList = getLineList($reader, 4)) {
		if(scalar(@lineList) == 4) {
			if((my $readLength = length($lineList[1])) == length($lineList[3])) {
				$readLengthCountHash{$readLength} += 1;
			} else {
				$error = 'different length';
				last;
			}
		} else {
			$error = 'insufficient line';
			last;
		}
	}
	close($reader);
	if($error eq '') {
		foreach my $readLength (sort {$b <=> $a} keys %readLengthCountHash) {
			print join("\t", $fastqFile, $readLength, $readLengthCountHash{$readLength}), "\n";
		}
	} else {
		print STDERR join("\t", $fastqFile, $error), "\n";
	}
}

sub getLineList {
	my ($reader, $number) = @_;
	my @lineList = ();
	foreach(1 .. $number) {
		if(defined(my $line = <$reader>)) {
			chomp($line);
			push(@lineList, $line);
		}
	}
	return @lineList;
}
