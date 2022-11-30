#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;

my ($antismashOutputDirectory) = @ARGV;
my @columnList = ('anchor', 'seq_id', 'start', 'end', 'tool', 'neighbouring_start', 'neighbouring_end', 'product', 'kind', 'prefix', 'height');
my %columnHash = map {$_ => 1} @columnList;
my $json = '';
{
	open(my $reader, "$antismashOutputDirectory/regions.js");
	while(my $line = <$reader>) {
		chomp($line);
		$line =~ s/^var (\S+) = /"$1": /;
		$line =~ s/;$/,/;
		$json .= $line;
	}
	close($reader);
	$json =~ s/,$//;
	$json = "{$json}";
}
$json = decode_json($json);
my @clusterList = ();
foreach my $record (@{$json->{'recordData'}}) {
	foreach my $region (@{$record->{'regions'}}) {
		foreach my $cluster (@{$region->{'clusters'}}) {
			$cluster->{'anchor'} = $region->{'anchor'};
			$cluster->{'seq_id'} = $record->{'seq_id'};
			push(@clusterList, $cluster);
		}
	}
}
foreach my $cluster (@clusterList) {
	foreach my $column (keys %$cluster) {
		next if($columnHash{$column});
		push(@columnList, $column);
		$columnHash{$column} = 1;
	}
}
%columnHash = ();
foreach my $cluster (@clusterList) {
	foreach my $column (keys %$cluster) {
		$columnHash{$column} = 1;
	}
}
@columnList = grep {$columnHash{$_}} @columnList;
print join("\t", @columnList), "\n";
foreach my $cluster (@clusterList) {
	print join("\t", map {defined($_) ? $_ : ''} @$cluster{@columnList}), "\n";
}
