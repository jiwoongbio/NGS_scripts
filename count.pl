#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $header = ''),
);
my ($inputFile, %countHash) = @ARGV;

open(my $reader, $inputFile);
if($header) {
	chomp(my $line = <$reader>);
	print join("\t", $line, 'count'), "\n";
}
while(my $line = <$reader>) {
	chomp($line);
	$countHash{$line} += 1;
}
close($reader);

print join("\t", $_, $countHash{$_}), "\n" foreach(sort keys %countHash);
