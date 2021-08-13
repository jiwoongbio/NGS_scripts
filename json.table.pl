#!/bin/env perl
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
$json = decode_json($json);
foreach my $key (@keyList) {
	$json = $json->{$key};
}
if(ref($json) eq 'HASH') {
	my @keyList = sort keys %$json;
	my %columnHash = ();
	foreach my $key (@keyList) {
		$columnHash{$_} = 1 foreach(keys %{$json->{$key}});
	}
	my @columnList = sort keys %columnHash;
	print join("\t", '', @columnList), "\n";
	foreach my $key (@keyList) {
		print join("\t", $key, map {defined($_) ? encode_json($_) : ''} map {$json->{$key}->{$_}} @columnList), "\n";
	}
}
