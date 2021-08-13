#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(max min);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   paf.minimum_match.pl [options] input.paf minimum_match

Options: -h       display this help message

EOF
}
my ($pafFile, $minimumMatch) = @ARGV;
my @pafMandatoryColumnList = ('query_name', 'query_length', 'query_start', 'query_end', 'strand', 'target_name', 'target_length', 'target_start', 'target_end', 'match', 'alignment_length', 'mapping_quality');
open(my $reader, ($pafFile =~ /\.gz$/) ? "gzip -dc $pafFile |" : $pafFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^@/) {
		print "$line\n";
	} else {
		my %tokenHash = ();
		(@tokenHash{@pafMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line, -1);
#		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		my @lengthList = @tokenHash{'query_length', 'target_length'};
		if($tokenHash{'strand'} eq '+') {
			push(@lengthList, max($tokenHash{'query_length'} - $tokenHash{'query_start'}, $tokenHash{'target_end'}));
			push(@lengthList, max($tokenHash{'query_end'}, $tokenHash{'target_length'} - $tokenHash{'target_start'}));
		}
		if($tokenHash{'strand'} eq '+') {
			push(@lengthList, max($tokenHash{'query_length'} - $tokenHash{'query_start'}, $tokenHash{'target_length'} - $tokenHash{'target_start'}));
			push(@lengthList, max($tokenHash{'query_end'}, $tokenHash{'target_end'}));
		}
		print "$line\n" if($tokenHash{'match'} / min(@lengthList) >= $minimumMatch);
	}
}
close($reader);
