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
	'v' => \(my $invert = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] read_name.txt input.fastq.gz | gzip > extract.fastq.gz

Options: -h       display this help message
         -v       invert

EOF
}
my ($readNameFile, $fastqFile) = @ARGV;
{
	my $pid = open2(my $reader, my $writer, "LC_ALL=C sort -t '\t' -k1,1");
	{
		open(my $reader, ($readNameFile =~ /\.gz$/ ? "gzip -dc $readNameFile |" : $readNameFile));
		while(my $line = <$reader>) {
			chomp($line);
			$line =~ s/\t/ /g;
			my $readName = $line;
			$readName =~ s/ .*$//;
			$readName =~ s/\/1$//;
			$readName =~ s/\/2$//;
			print $writer "$readName\n";
		}
		close($reader);
		close($writer);
	}
	sub getReadName {
		while(my $line = <$reader>) {
			chomp($line);
			my $readName = $line;
			return $readName;
		}
		return '';
	}
	sub closeReadNameReader {
		close($reader);
		waitpid($pid, 0);
	}
}
{
	my $pid = open2(my $reader, my $writer, "LC_ALL=C sort -t '\t' -k1,1");
	{
		open(my $reader, ($fastqFile =~ /\.gz$/ ? "gzip -dc $fastqFile |" : $fastqFile));
		while(my @lineList = getLineList($reader, 4)) {
			if(scalar(@lineList) == 4) {
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
	my $readName = getReadName();
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		$readName = getReadName() while($readName ne '' && $readName lt $tokenList[0]);
		if($invert eq '') {
			if($readName ne '' && $readName eq $tokenList[0]) {
				print "$_\n" foreach(@tokenList[1 .. 4]);
			}
		} else {
			unless($readName ne '' && $readName eq $tokenList[0]) {
				print "$_\n" foreach(@tokenList[1 .. 4]);
			}
		}
	}
	close($reader);
	waitpid($pid, 0);
}
closeReadNameReader();

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
