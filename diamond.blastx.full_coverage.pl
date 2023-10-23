#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my @samMandatoryColumnList = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');
my ($samFile) = @ARGV;
open(my $reader, ($samFile =~ /\.gz$/ ? "gzip -dc $samFile |" : $samFile));
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^@/) {
		print "$line\n";
	} else {
		my %tokenHash = ();
		(@tokenHash{@samMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		$tokenHash{'qname'} =~ s/_([0-9]+)$//;
		my $readLength = $1;
		if($tokenHash{'MD:Z'} =~ /^[0-9]+$/ || ($tokenHash{'pos'} == 1 && $tokenHash{'MD:Z'} =~ /^[A-Z][0-9]+$/)) {
			if($tokenHash{'ZF:i'} < 0) {
				if($tokenHash{'ZL:i'} == length($tokenHash{'seq'})) {
					print "$line\n";
				} elsif($tokenHash{'ZS:i'} == $readLength) {
					print "$line\n";
				} elsif($tokenHash{'ZS:i'} - length($tokenHash{'seq'}) * 3 < 3) {
					print "$line\n";
				}
			}
		}
	}
}
close($reader);
