#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'r' => \(my $ignoreReferenceMatch = ''),
);
my (@samFileList) = @ARGV;
my %strandReadCountHash = ();
foreach my $samFile (@samFileList) {
	open(my $reader, "samtools view -F 2316 $samFile |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		if($tokenHash{'flag'} & 1) {
			next if($ignoreReferenceMatch eq '' && $tokenHash{'rnext'} ne '=');
			$strandReadCountHash{'reverse'} += 1 if(($tokenHash{'flag'} & 240) == 80 || ($tokenHash{'flag'} & 240) == 160);
			$strandReadCountHash{'forward'} += 1 if(($tokenHash{'flag'} & 240) == 96 || ($tokenHash{'flag'} & 240) == 144);
		} else {
			$strandReadCountHash{'reverse'} += 1 if(($tokenHash{'flag'} & 240) == 16);
			$strandReadCountHash{'forward'} += 1 if(($tokenHash{'flag'} & 240) == 0);
		}
	}
}
foreach my $strand ('reverse', 'forward') {
	my $readCount = $strandReadCountHash{$strand};
	$readCount = 0 unless(defined($readCount));
	print join("\t", $strand, $readCount), "\n";
}
