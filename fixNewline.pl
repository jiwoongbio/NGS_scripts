#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

foreach my $file (@ARGV) {
	my $fixed = '';
	open(my $reader, $file);
	open(my $writer, "> $file.newlineFixed");
	while(my $line = <$reader>) {
		$fixed = 1 if($line !~ /\n$/ && $line =~ s/$/\n/);
		$fixed = 1 if($line =~ s/\r\n?/\n/g);
		print $writer $line;
	}
	close($reader);
	close($writer);
	if($fixed) {
		system("mv $file.newlineFixed $file");
	} else {
		system("rm $file.newlineFixed");
	}
}
