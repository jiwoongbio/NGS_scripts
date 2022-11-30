#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min max);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'S=s' => \(my $stranded = ''),
);
my ($bamFile, $chromosome, $start, $end) = @ARGV;
open(my $reader, "samtools view -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$end-$end |");
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
	$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
	next if($tokenHash{'pos'} > $start);
	if($stranded eq 'f' || $stranded eq 'forward') {
		next unless(($tokenHash{'flag'} & 253) == 97 || ($tokenHash{'flag'} & 253) == 145 || ($tokenHash{'flag'} & 253) == 0);
	}
	if($stranded eq 'r' || $stranded eq 'reverse') {
		next unless(($tokenHash{'flag'} & 253) == 81 || ($tokenHash{'flag'} & 253) == 161 || ($tokenHash{'flag'} & 253) == 16);
	}
	my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
	my $index = max(0, map {$_ + 1} grep {$positionList[$_] ne '' && $positionList[$_] < $start} 0 .. $#positionList);
	my $length = min(scalar(@positionList), grep {$positionList[$_] ne '' && $positionList[$_] > $end} 0 .. $#positionList) - $index;
	my $name = $tokenHash{'qname'};
	my $number = ($tokenHash{'flag'} & 192) / 64;
	$name = "$name/$number" if($number);
	print join("\t", $name, substr($tokenHash{'seq'}, $index, $length)), "\n";
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
