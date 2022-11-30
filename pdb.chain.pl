#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($pdbFile) = @ARGV;
my %chainHash = ();
open(my $reader, ($pdbFile =~ /\.gz$/ ? "gzip -dc $pdbFile |" : $pdbFile));
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^ATOM  /) {
		my $chain = substr($line, 21, 1);
		$chainHash{$chain} = 1;
	}
}
close($reader);

foreach my $chain (sort keys %chainHash) {
	print "$chain\n";
}
