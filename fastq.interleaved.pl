#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] input.1.fastq.gz input.2.fastq.gz | gzip > interleaved.fastq.gz

Options: -h       display this help message

EOF
}
my ($fastqFile1, $fastqFile2) = @ARGV;
open(my $reader1, ($fastqFile1 =~ /\.gz$/ ? "gzip -dc $fastqFile1 |" : $fastqFile1));
open(my $reader2, ($fastqFile2 =~ /\.gz$/ ? "gzip -dc $fastqFile2 |" : $fastqFile2));
my @lineList1 = getLineList($reader1, 4);
my @lineList2 = getLineList($reader2, 4);
while(@lineList1 || @lineList2) {
	die "wrong format fastq 1\n" unless(scalar(@lineList1) == 4 && length($lineList1[1]) == length($lineList1[3]) && $lineList1[0] =~ /^\@/ && $lineList1[2] =~ /^\+/);
	die "wrong format fastq 2\n" unless(scalar(@lineList2) == 4 && length($lineList2[1]) == length($lineList2[3]) && $lineList2[0] =~ /^\@/ && $lineList2[2] =~ /^\+/);
	$lineList1[0] =~ s/\s.*$//;
	$lineList2[0] =~ s/\s.*$//;
	$lineList1[0] =~ s/\/1$//;
	$lineList2[0] =~ s/\/2$//;
	die "different sequence name\n" unless($lineList1[0] eq $lineList2[0]);
	$lineList1[0] .= '/1';
	$lineList2[0] .= '/2';
	$lineList1[2] = '+';
	$lineList2[2] = '+';
	print "$_\n" foreach(@lineList1, @lineList2);
	@lineList1 = getLineList($reader1, 4);
	@lineList2 = getLineList($reader2, 4);
}
close($reader1);
close($reader2);

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
