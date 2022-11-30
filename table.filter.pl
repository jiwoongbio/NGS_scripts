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
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] table.txt column operator value > table.filtered.txt

Options: -h       display this help message

EOF
}
my ($tableFile, $column, $operator, $value) = @ARGV;
open(my $reader, ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile |" : $tableFile));
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line, -1);
print join("\t", @columnList), "\n";
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my %tokenHash = ();
	@tokenHash{@columnList} = @tokenList;
	if($operator eq '=~' || $operator eq '!~') {
		if(eval(sprintf('$tokenHash{$column} %s /$value/', $operator))) {
			print join("\t", @tokenList), "\n";
		}
	} elsif($operator =~ /^[!=<>]+$/) {
		if(eval(sprintf('$tokenHash{$column} %s $value', $operator))) {
			print join("\t", @tokenList), "\n";
		}
	} else {
		if(eval(sprintf('$tokenHash{$column} %s "$value"', $operator))) {
			print join("\t", @tokenList), "\n";
		}
	}
}
close($reader);
