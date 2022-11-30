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
	'read' => \(my $printReadVariant = ''),
);
my (@samFileList) = @ARGV;
foreach my $samFile (@samFileList) {
	open(my $reader, "samtools view -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $samFile |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line, -1);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		my ($read, $number) = ($tokenHash{'qname'}, ($tokenHash{'flag'} & 192) / 64);
		next if($tokenHash{'flag'} & 4);
		foreach(getVariantList(@tokenHash{'rname', 'pos', 'cigar', 'MD:Z', 'seq'})) {
			my ($chromosome, $position, $refBase, $altBase, $index) = @$_;
			if($printReadVariant eq '') {
				print join("\t", $chromosome, $position, $refBase, $altBase, $read, $number), "\n";
			} elsif(($tokenHash{'flag'} & 16) == 0) {
				print join("\t", (($number ? "$read/$number" : $read), $index + 1, $altBase, $refBase)), "\n";
			} else {
				print join("\t", (($number ? "$read/$number" : $read), length($tokenHash{'seq'}) - ($index + length($altBase)) + 1, getReverseComplementarySequence($altBase), getReverseComplementarySequence($refBase))), "\n";
			}
		}
	}
	close($reader);
}

sub getVariantList {
	my ($chromosome, $position, $cigar, $md, $sequence) = @_;
	my @variantList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			if($md =~ s/^([0-9]+)//) {
				if($length > $1) {
					$cigar = join('', $length - $1, 'M', $cigar);
					$length = $1;
				} elsif($1 > $length) {
					$md = join('', $1 - $length, $md);
				}
				$index += $length;
				$position += $length;
			} elsif($md =~ s/^([A-Z0]+)//) {
				(my $refBase = $1) =~ s/0//g;
				if($length > length($refBase)) {
					$cigar = join('', $length - length($refBase), 'M', $cigar);
					$length = length($refBase);
				} elsif(length($refBase) > $length) {
					$md = join('', substr($refBase, $length), $md);
					$refBase = substr($refBase, 0, $length);
				}
				my $altBase = substr($sequence, $index, $length);
				push(@variantList, [$chromosome, $position, $refBase, $altBase, $index]);
				$index += $length;
				$position += $length;
			} else {
				die "Failed to parse CIGAR/MD";
			}
		} elsif($operation eq 'I') {
				my ($refBase, $altBase) = ('', substr($sequence, $index, $length));
				push(@variantList, [$chromosome, $position, $refBase, $altBase, $index]);
				$index += $length;
		} elsif($operation eq 'D') {
			if($md =~ s/^\^([A-Z]+)0*// && length($1) == $length) {
				my ($refBase, $altBase) = ($1, '');
				push(@variantList, [$chromosome, $position, $refBase, $altBase, $index]);
				$position += $length;
			} else {
				die "Failed to parse CIGAR/MD";
			}
		} elsif($operation eq 'N') {
				$position += $length;
		} elsif($operation eq 'S') {
				$index += $length;
		} elsif($operation eq 'H') {
		} else {
				die "Failed to parse CIGAR/MD";
		}
	}
	die "Failed to parse CIGAR/MD" unless($cigar eq '' && $md eq '' && $index == length($sequence));
	return @variantList;
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
