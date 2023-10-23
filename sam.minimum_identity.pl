#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum min max all);
use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
	'd' => \(my $isSamFromDiamond = ''),
	'U' => \(my $ignoreUnmappedRead = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] input.sam.gz minimum_identity | gzip > filtered.sam.gz

Options: -h       display this help message
         -d       is sam from diamond
         -U       ignore unmapped read

EOF
}
my ($samFile, $minimumIdentity) = @ARGV;
open(my $reader, ($samFile =~ /\.gz$/) ? "gzip -dc $samFile |" : $samFile);
my %sequenceLengthHash = ();
my @tokenListList = ();
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	if($line =~ /^@/) {
		if($tokenList[0] eq '@SQ') {
			my %tokenHash = map {$_->[0] => $_->[1]} map {[split(/:/, $_, 2)]} @tokenList;
			$sequenceLengthHash{$tokenHash{'SN'}} = $tokenHash{'LN'};
		}
		print "$line\n";
	} else {
		if(@tokenListList && $tokenListList[0]->[0] ne $tokenList[0]) {
			if(all {$_ >= $minimumIdentity} getMaximumIdentityList(@tokenListList)) {
				print join("\t", @$_), "\n" foreach(@tokenListList);
			}
			@tokenListList = ();
		}
		push(@tokenListList, \@tokenList);
	}
}
if(@tokenListList) {
	if(all {$_ >= $minimumIdentity} getMaximumIdentityList(@tokenListList)) {
		print join("\t", @$_), "\n" foreach(@tokenListList);
	}
	@tokenListList = ();
}
close($reader);

sub getMaximumIdentityList {
	my @tokenListList = @_;
	my %readNumberIdentityListHash = ();
	push(@{$readNumberIdentityListHash{($_->[1] & 192) / 64}}, getIdentity(@$_)) foreach(grep {$ignoreUnmappedRead eq '' || ($_->[1] & 4) == 0} @tokenListList);
	return map {max(@{$readNumberIdentityListHash{$_}})} sort {$a <=> $b} keys %readNumberIdentityListHash;
}

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
