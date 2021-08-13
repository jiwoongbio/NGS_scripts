#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($fastaFile, $lineLength) = @ARGV;
my $sequence = '';
open(my $reader, $fastaFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^>/) {
		printSequence($sequence) if($sequence ne '');
		print "$line\n";
		$sequence = '';
	} else {
		$sequence .= $line;
	}
}
close($reader);
printSequence($sequence) if($sequence ne '');

sub printSequence {
	my ($sequence) = @_;
	for(my $index = 0; $index < length($sequence); $index += $lineLength) {
		print substr($sequence, $index, $lineLength), "\n";
	}
}
