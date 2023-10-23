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
	'n=s' => \(my $nameColumnIndexes = 0),
	'v=s' => \(my $valueColumnIndexes = 1),
	's=i' => \(my $skipLines = 0),
	'd=s' => \(my $defaultValue = 0),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   $file [options] column=files [...]

Options: -h       display this help message
         -n STR   name column indexes separated by comma [$nameColumnIndexes]
         -v STR   value column indexes separated by comma [$valueColumnIndexes]
         -s INT   skip first INT lines [$skipLines]
         -d STR   degault value [$defaultValue]

EOF
}
my @nameColumnIndexList = eval($nameColumnIndexes);
my @valueColumnIndexList = eval($valueColumnIndexes);

my @columnFilesList = @ARGV;
@columnFilesList = map {[split(/=/, $_, 2)]} @columnFilesList;

my %nameValueListHash = ();
foreach my $index (0 .. $#columnFilesList) {
	my ($column, $files) = @{$columnFilesList[$index]}[0, -1];
	foreach my $file (split(/,/, $files)) {
		open(my $reader, $file);
		foreach(1 .. $skipLines) {
			my $line = <$reader>;
		}
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line, -1);
			my $name = join("\t", @tokenList[@nameColumnIndexList]);
			if($nameColumnIndexes eq $valueColumnIndexes) {
				$nameValueListHash{$name}->[$index] = 1;
			} elsif(scalar(@valueColumnIndexList) == 1) {
				if(defined($nameValueListHash{$name}->[$index])) {
					$nameValueListHash{$name}->[$index] += $tokenList[$valueColumnIndexList[0]];
				} else {
					$nameValueListHash{$name}->[$index] = $tokenList[$valueColumnIndexList[0]];
				}
			} else {
				$nameValueListHash{$name}->[$index] = join(',', @tokenList[@valueColumnIndexList]);
			}
		}
		close($reader);
	}
}

print join("\t", ('') x scalar(@nameColumnIndexList), map {$_->[0]} @columnFilesList), "\n";
foreach my $name (sort keys %nameValueListHash) {
	print join("\t", $name, map {defined($_) ? $_ : $defaultValue} @{$nameValueListHash{$name}}[0 .. $#columnFilesList]), "\n";
}
