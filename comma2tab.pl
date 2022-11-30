#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($commaFile) = @ARGV;
open(my $reader, $commaFile);
chomp(my $line = <$reader>);
my @lastTokenList = split(/,/, $line, -1);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/,/, $line, -1);
	if(@lastTokenList && $lastTokenList[-1] =~ /^"/ && $lastTokenList[-1] !~ /"$/) { # start double quote, but not end double quote
		if(@tokenList) {
			$lastTokenList[-1] = join(',', $lastTokenList[-1], $tokenList[0]); # merge last token and first token
			push(@lastTokenList, @tokenList[1 .. $#tokenList]); # merge with next line
		}
	} else {
		print join("\t", remergeDoubleQuotes(@lastTokenList)), "\n";
		@lastTokenList = @tokenList;
	}
}
print join("\t", remergeDoubleQuotes(@lastTokenList)), "\n";
close($reader);

sub remergeDoubleQuotes {
	my @tokenList = @_;
	for(my $index = 0; $index + 1 < scalar(@tokenList);) {
		if($tokenList[$index] =~ /^"/ && $tokenList[$index] !~ /"$/) { # start double quote, but not end double quote
			splice(@tokenList, $index, 2, join(',', @tokenList[$index, $index + 1])); # merge with next token
		} else {
			$index += 1;
		}
	}
	@tokenList = map {$_ =~ /^"(.*)"$/ ? $1 : $_} @tokenList; # remove double quotes
	return @tokenList;
}
