#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
	'H' => \(my $hasHeader = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] input.txt [key1 count1 ...] > key.count.txt

Options: -h       display this help message
         -H       has a header line

EOF
}
my ($inputFile, %countHash) = @ARGV;

open(my $reader, $inputFile);
if($hasHeader) {
	chomp(my $line = <$reader>);
	print join("\t", $line, 'count'), "\n";
}
while(my $line = <$reader>) {
	chomp($line);
	$countHash{$line} += 1;
}
close($reader);

print join("\t", $_, $countHash{$_}), "\n" foreach(sort keys %countHash);
