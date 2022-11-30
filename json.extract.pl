#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;

my ($jsonFile, @keyList) = @ARGV;
my $json = '';
{
	open(my $reader, $jsonFile);
	$json .= $_ while(<$reader>);
	close($reader);
}
my @jsonList = (decode_json($json));
foreach my $key (@keyList) {
	@jsonList = map {ref($_) eq 'ARRAY' ? @$_ : $_} @jsonList while(grep {ref($_) eq 'ARRAY'} @jsonList);
	@jsonList = map {defined($_) ? $_ : ()} map {$_->{$key}} @jsonList;
}
@jsonList = map {ref($_) eq 'ARRAY' ? @$_ : $_} @jsonList while(grep {ref($_) eq 'ARRAY'} @jsonList);
print "$_\n" foreach(@jsonList);
