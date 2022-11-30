#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'o' => \(my $objectRow = ''),
	'n' => \(my $numberValue = ''),
);
my ($tableFile) = @ARGV;
my @rowList = ();
open(my $reader, ($tableFile =~ /\.gz$/ ? "gzip -dc $tableFile |" : $tableFile));
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line, -1);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	if($objectRow eq '') {
		if($numberValue eq '') {
			push(@rowList, sprintf('[%s]', join(',', sprintf('"%s"', $tokenList[0]), map {sprintf('"%s"', $tokenList[$_])} 1 .. $#tokenList)));
		} else {
			push(@rowList, sprintf('[%s]', join(',', sprintf('"%s"', $tokenList[0]), map {sprintf('%s', $tokenList[$_])} 1 .. $#tokenList)));
		}
	} else {
		if($numberValue eq '') {
			push(@rowList, sprintf('{%s}', join(',', sprintf('"%s":"%s"', $columnList[0], $tokenList[0]), map {$tokenList[$_] ? sprintf('"%s":"%s"', $columnList[$_], $tokenList[$_]) : ()} 1 .. $#tokenList)));
		} else {
			push(@rowList, sprintf('{%s}', join(',', sprintf('"%s":"%s"', $columnList[0], $tokenList[0]), map {$tokenList[$_] ? sprintf('"%s":%s', $columnList[$_], $tokenList[$_]) : ()} 1 .. $#tokenList)));
		}
	}
}
close($reader);
printf('{"columns":[%s],"rows":[%s]}', join(',', map {sprintf('"%s"', $_)} @columnList), join(',', @rowList));
