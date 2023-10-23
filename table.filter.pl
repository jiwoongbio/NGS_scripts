#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(all);
use Getopt::Long qw(:config no_ignore_case);

chomp(my $file = `readlink -f $0`);
my $directory = $1 if($file =~ s/(^.*\/)//);

GetOptions(
	'h' => \(my $help = ''),
	'b=s' => \(my $booleanOperator = 'and'),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] table.txt column operator value [...] > table.filtered.txt

Options: -h       display this help message
         -b STR   boolean operator [$booleanOperator]

EOF
}
my ($tableFile) = @ARGV;
my @columnOperatorValueList = ();
for(my $index = 1; $index < scalar(@ARGV); $index += 3) {
	my ($column, $operator, $value) = @ARGV[$index .. $index + 2];
	push(@columnOperatorValueList, [$column, $operator, $value]);
}
open(my $reader, ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile |" : $tableFile));
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line, -1);
print join("\t", @columnList), "\n";
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my %tokenHash = ();
	@tokenHash{@columnList} = @tokenList;
	my @booleanList = ();
	foreach(@columnOperatorValueList) {
		my ($column, $operator, $value) = @$_;
		if($operator eq '=~' || $operator eq '!~') {
			push(@booleanList, eval(sprintf('$tokenHash{$column} %s /$value/', $operator)));
		} elsif($operator =~ /^[!=<>]+$/) {
			push(@booleanList, eval(sprintf('$tokenHash{$column} %s $value', $operator)));
		} else {
			push(@booleanList, eval(sprintf('$tokenHash{$column} %s "$value"', $operator)));
		}
	}
	if($booleanOperator eq 'and') {
		print join("\t", @tokenList), "\n" if(all {$_} @booleanList[0 .. $#columnOperatorValueList]);
	}
	if($booleanOperator eq 'or') {
		print join("\t", @tokenList), "\n" if(grep {$_} @booleanList[0 .. $#columnOperatorValueList]);
	}
}
close($reader);
