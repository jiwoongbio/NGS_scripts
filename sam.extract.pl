#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'E' => \(my $noEmpty = ''),
);
my @samMandatoryColumnList = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');
my ($samFile, @columnList) = @ARGV;
open(my $reader, ($samFile =~ /\.gz$/ ? "gzip -dc $samFile |" : $samFile));
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^@/);
	my %tokenHash = ();
	(@tokenHash{@samMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line);
	$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
	next if($tokenHash{'mapq'} < $minimumMappingQuality);
	next if((int($tokenHash{'flag'}) & $includeFlag) != $includeFlag);
	next if((int($tokenHash{'flag'}) & $excludeFlag) != 0);
	my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
	@positionList = grep {$_ ne ''} @positionList;
	@tokenHash{'chromosome', 'start', 'end', 'strand'} = ($tokenHash{'rname'}, @positionList[0, -1], (($tokenHash{'flag'} & 16) == 0 ? '+' : '-'));
	my @tokenList = @tokenHash{@columnList};
	next if($noEmpty && grep {!defined($_)} @tokenList);
	print join("\t", map {defined($_) ? $_ : ''} @tokenList), "\n";
}
close($reader);

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		} elsif($operation eq 'D') {
			$position += $length;
		} elsif($operation eq 'N') {
			$position += $length;
		} elsif($operation eq 'S') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		}
	}
	return @positionList;
}
