#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);
use IPC::Open2;

my @samMandatoryColumnList = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');
my ($samFile) = @ARGV;

my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2n -k3,3n");
{
	open(my $reader, ($samFile =~ /\.gz$/ ? "gzip -dc $samFile |" : $samFile));
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^@/) {
		} else {
			my %tokenHash = ();
			(@tokenHash{@samMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line);
			$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
			$tokenHash{'qname'} =~ s/_([0-9]+)$//;
			my $readLength = $1;
			print $writer join("\t", $tokenHash{'rname'}, $tokenHash{'pos'}, $tokenHash{'pos'} + length($tokenHash{'seq'}) - 1, $tokenHash{'ZL:i'}), "\n";
		}
	}
	close($reader);
}
close($writer);
{
	my ($name, $length, $count, @depthList) = ('', 0);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{'name', 'start', 'end', 'length'} = split(/\t/, $line, -1);
		if($tokenHash{'name'} ne $name) {
			if($length > 0) {
				print join("\t", $name, $count, sum(0, grep {defined} @depthList) / $length, scalar(grep {defined} @depthList) / $length), "\n";
			}
			($name, $length, $count, @depthList) = (@tokenHash{'name', 'length'}, 0);
		}
		$count += 1;
		$depthList[$_] += 1 foreach(map {$_ - 1} $tokenHash{'start'} .. $tokenHash{'end'});
	}
	if($length > 0) {
		print join("\t", $name, $count, sum(0, grep {defined} @depthList) / $length, scalar(grep {defined} @depthList) / $length), "\n";
	}
}
close($reader);
waitpid($pid, 0);
