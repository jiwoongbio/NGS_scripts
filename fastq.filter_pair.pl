#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] input.1.fastq.gz input.2.fastq.gz output.1.fastq.gz output.2.fastq.gz

Options: -h       display this help message

EOF
}
my ($inputFastqFile1, $inputFastqFile2, $outputFastqFile1, $outputFastqFile2) = @ARGV;
{
	my $pid = open2(my $reader, my $writer, "LC_ALL=C sort -t '\t' -k1,1 -k2 | uniq");
	{
		open(my $reader, ($inputFastqFile1 =~ /\.gz$/ ? "gzip -dc $inputFastqFile1 |" : $inputFastqFile1));
		while(my @lineList = getLineList($reader, 4)) {
			if(scalar(@lineList) == 4 && length($lineList[1]) == length($lineList[3]) && $lineList[0] =~ /^\@/ && $lineList[2] =~ /^\+/) {
				s/\t/ /g foreach(@lineList);
				(my $readName = $lineList[0]) =~ s/^\@//;
				$readName =~ s/ .*$//;
				$readName =~ s/\/1$//;
				$readName =~ s/\/2$//;
				print $writer join("\t", $readName, @lineList), "\n";
			}
		}
		close($reader);
		close($writer);
	}
	sub readNameTokenList1 {
		while(my $line = <$reader>) {
			chomp($line);
			my @readNameTokenList = split(/\t/, $line);
			return @readNameTokenList;
		}
		return ();
	}
	sub closeReader1 {
		close($reader);
		waitpid($pid, 0);
	}
}
{
	my $pid = open2(my $reader, my $writer, "LC_ALL=C sort -t '\t' -k1,1 -k2 | uniq");
	{
		open(my $reader, ($inputFastqFile2 =~ /\.gz$/ ? "gzip -dc $inputFastqFile2 |" : $inputFastqFile2));
		while(my @lineList = getLineList($reader, 4)) {
			if(scalar(@lineList) == 4 && length($lineList[1]) == length($lineList[3]) && $lineList[0] =~ /^\@/ && $lineList[2] =~ /^\+/) {
				s/\t/ /g foreach(@lineList);
				(my $readName = $lineList[0]) =~ s/^\@//;
				$readName =~ s/ .*$//;
				$readName =~ s/\/1$//;
				$readName =~ s/\/2$//;
				print $writer join("\t", $readName, @lineList), "\n";
			}
		}
		close($reader);
		close($writer);
	}
	sub readNameTokenList2 {
		while(my $line = <$reader>) {
			chomp($line);
			my @readNameTokenList = split(/\t/, $line);
			return @readNameTokenList;
		}
		return ();
	}
	sub closeReader2 {
		close($reader);
		waitpid($pid, 0);
	}
}
{
	open(my $writer, ($outputFastqFile1 =~ /\.gz$/ ? "| gzip > $outputFastqFile1" : "> $outputFastqFile1"));
	sub writeTokenList1 {
		print $writer "$_\n" foreach(@_);
	}
	sub closeWriter1 {
		close($writer);
	}
}
{
	open(my $writer, ($outputFastqFile2 =~ /\.gz$/ ? "| gzip > $outputFastqFile2" : "> $outputFastqFile2"));
	sub writeTokenList2 {
		print $writer "$_\n" foreach(@_);
	}
	sub closeWriter2 {
		close($writer);
	}
}
my $count = 0;
my @readNameTokenList1 = readNameTokenList1();
my @readNameTokenList2 = readNameTokenList2();
while(@readNameTokenList1 && @readNameTokenList2) {
	if($readNameTokenList1[0] eq $readNameTokenList2[0]) {
		writeTokenList1(@readNameTokenList1[1 .. 4]);
		writeTokenList2(@readNameTokenList2[1 .. 4]);
		@readNameTokenList1 = readNameTokenList1();
		@readNameTokenList2 = readNameTokenList2();
		$count += 1;
	} elsif($readNameTokenList1[0] lt $readNameTokenList2[0]) {
		@readNameTokenList1 = readNameTokenList1();
	} elsif($readNameTokenList1[0] gt $readNameTokenList2[0]) {
		@readNameTokenList2 = readNameTokenList2();
	}
}
closeReader1();
closeReader2();
closeWriter1();
closeWriter2();

print "$count\n";

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
