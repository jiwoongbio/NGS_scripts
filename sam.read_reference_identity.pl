#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum min all);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'd' => \(my $isSamFromDiamond = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   sam.read_identity.pl [options] input.sam

Options: -h       display this help message
         -d       is sam from diamond

EOF
}
my ($samFile) = @ARGV;
open(my $reader, ($samFile =~ /\.gz$/) ? "gzip -dc $samFile |" : $samFile);
my %sequenceLengthHash = ();
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	if($line =~ /^@/) {
		if($tokenList[0] eq '@SQ') {
			my %tokenHash = map {$_->[0] => $_->[1]} map {[split(/:/, $_, 2)]} @tokenList;
			$sequenceLengthHash{$tokenHash{'SN'}} = $tokenHash{'LN'};
		}
	} else {
		print join("\t", $tokenList[0], ($tokenList[1] & 192) / 64, $tokenList[2], getIdentity(@tokenList)), "\n";
	}
}
close($reader);

sub getIdentity {
	my %tokenHash = ();
	(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = @_;
	$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
	return 0 if($tokenHash{'flag'} & 4);
	my $sequenceLength = $sequenceLengthHash{$tokenHash{'rname'}};
	$sequenceLength = $tokenHash{'ZL:i'} if($isSamFromDiamond);
	my $identity = sum(0, ($tokenHash{'MD:Z'} =~ /([0-9]+)/g)) / min(
		sum(length($tokenHash{'seq'}), ($tokenHash{'cigar'} =~ /([0-9]+)D/g)), # reference - alignment - reference
		sum(($tokenHash{'pos'} - 1), ($tokenHash{'cigar'} =~ /([0-9]+)[MDN]/g), ($tokenHash{'cigar'} =~ /([0-9]+)S$/), ($tokenHash{'cigar'} =~ /([0-9]+)I/g)), # read - alignment - reference
		sum($sequenceLength - ($tokenHash{'pos'} - 1), ($tokenHash{'cigar'} =~ /^([0-9]+)S/), ($tokenHash{'cigar'} =~ /([0-9]+)I/g)), # reference - alignment - read
		sum($sequenceLength, ($tokenHash{'cigar'} =~ /([0-9]+)I/g)), # read - alignment - read
	);
	return $identity;
}
