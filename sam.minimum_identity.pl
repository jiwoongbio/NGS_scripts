#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum min);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   sam.minimum_identity.pl [options] input.sam minimum_identity

Options: -h       display this help message

EOF
}
my ($samFile, $minimumIdentity) = @ARGV;
open(my $reader, ($samFile =~ /\.gz$/) ? "gzip -dc $samFile |" : $samFile);
my %sequenceLengthHash = ();
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^@/) {
		my @tokenList = split(/\t/, $line, -1);
		if($tokenList[0] eq '@SQ') {
			my %tokenHash = map {$_->[0] => $_->[1]} map {[split(/:/, $_, 2)]} @tokenList;
			$sequenceLengthHash{$tokenHash{'SN'}} = $tokenHash{'LN'};
		}
		print "$line\n";
	} else {
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		next if($tokenHash{'flag'} & 4);
		my $identity = sum(0, ($tokenHash{'MD:Z'} =~ /([0-9]+)/g)) / min(
			sum(length($tokenHash{'seq'}), ($tokenHash{'cigar'} =~ /([0-9]+)D/g)), # reference - alignment - reference
			sum(($tokenHash{'pos'} - 1), ($tokenHash{'cigar'} =~ /([0-9]+)[MDN]/g), ($tokenHash{'cigar'} =~ /([0-9]+)S$/), ($tokenHash{'cigar'} =~ /([0-9]+)I/g)), # read - alignment - reference
			sum($sequenceLengthHash{$tokenHash{'rname'}} - ($tokenHash{'pos'} - 1), ($tokenHash{'cigar'} =~ /^([0-9]+)S/), ($tokenHash{'cigar'} =~ /([0-9]+)I/g)), # reference - alignment - read
			sum($sequenceLengthHash{$tokenHash{'rname'}}, ($tokenHash{'cigar'} =~ /([0-9]+)I/g)), # read - alignment - read
		);
		print "$line\n" if($identity >= $minimumIdentity);
	}
}
close($reader);
