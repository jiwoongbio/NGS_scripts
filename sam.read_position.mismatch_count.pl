#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
);
my (@samFileList) = @ARGV;
my $readCount = 0;
my %readPositionMismatchCountHash = ();
foreach my $samFile (@samFileList) {
	open(my $reader, "samtools view -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $samFile |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		next if($tokenHash{'flag'} & 4);
		next unless(defined($tokenHash{'MD:Z'}));
		$readCount += 1;
		$readPositionMismatchCountHash{$_} += 1 foreach(getMismatchReadPositionList(%tokenHash));
	}
	close($reader);
}
foreach my $readPosition (sort {$a <=> $b} grep {$_ > 0} keys %readPositionMismatchCountHash) {
	print join("\t", $readPosition, $readPositionMismatchCountHash{$readPosition}, $readPositionMismatchCountHash{$readPosition} / $readCount), "\n";
}
foreach my $readPosition (sort {$a <=> $b} grep {$_ < 0} keys %readPositionMismatchCountHash) {
	print join("\t", $readPosition, $readPositionMismatchCountHash{$readPosition}, $readPositionMismatchCountHash{$readPosition} / $readCount), "\n";
}

sub getMismatchReadPositionList {
	my (%tokenHash) = @_;
	my @indexLengthList = ();
	my ($cigar, $md, $position, $index) = (@tokenHash{'cigar', 'MD:Z', 'pos'}, 0);
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			if($md =~ s/^([0-9]+)// && $1 > 0) {
				if($length > $1) {
					$cigar = join('', $length - $1, 'M', $cigar);
					$length = $1;
				} elsif($1 > $length) {
					$md = join('', $1 - $length, $md);
				}
				$index += $length;
				$position += $length;
			} elsif($md =~ s/^([A-Z]+)0*//) {
				if($length > length($1)) {
					$cigar = join('', $length - length($1), 'M', $cigar);
					$length = length($1);
				} elsif(length($1) > $length) {
					$md = join('', substr($1, $length), $md);
				}
				push(@indexLengthList, [$index, $length]);
				$index += $length;
				$position += $length;
			} else {
				die "Failed to parse CIGAR/MD";
			}
		} elsif($operation eq 'I') {
				push(@indexLengthList, [$index, $length]);
				$index += $length;
		} elsif($operation eq 'D') {
			if($md =~ s/^\^([A-Z]+)0*// && length($1) == $length) {
				push(@indexLengthList, [$index, 0]);
				$position += $length;
			} else {
				die "Failed to parse CIGAR/MD";
			}
		} elsif($operation eq 'N') {
				$position += $length;
		} elsif($operation eq 'S') {
				push(@indexLengthList, [$index, $length]);
				$index += $length;
		} elsif($operation eq 'H') {
		} else {
				die "Failed to parse CIGAR/MD";
		}
	}
	die "Failed to parse CIGAR/MD" unless($cigar eq '' && $md eq '');
	my $length = length($tokenHash{'seq'});
	my @positionList = ();
	@indexLengthList = map {[$length - $_->[0] - $_->[1], $_->[1]]} reverse(@indexLengthList) if($tokenHash{'flag'} & 16);
	push(@positionList, map {"+$_"} map {$_->[1] == 0 ? ($_->[0] + 1) : ($_->[0] + 1 .. $_->[0] + $_->[1])} @indexLengthList);
	@indexLengthList = map {[$length - $_->[0] - $_->[1], $_->[1]]} reverse(@indexLengthList);
	push(@positionList, map {"-$_"} map {$_->[1] == 0 ? ($_->[0] + 1) : ($_->[0] + 1 .. $_->[0] + $_->[1])} @indexLengthList);
	return @positionList;
}
