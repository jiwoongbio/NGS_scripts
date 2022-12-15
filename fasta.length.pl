#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($fastaFile) = @ARGV;
my ($sequenceName, $sequenceLength) = ('', 0);
open(my $reader, $fastaFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^>(\S*)/) {
		print join("\t", $sequenceName, $sequenceLength), "\n" if($sequenceLength > 0);
		($sequenceName, $sequenceLength) = ($1, 0);
	} else {
		$sequenceLength += length($line);
	}
}
close($reader);
print join("\t", $sequenceName, $sequenceLength), "\n" if($sequenceLength > 0);
