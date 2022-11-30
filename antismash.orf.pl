#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;

my ($antismashOutputDirectory) = @ARGV;
my @columnList = ('anchor', 'seq_id', 'start', 'end', 'strand', 'locus_tag', 'type', 'description');
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
my @orfList = ();
foreach my $record (@{$json->{'recordData'}}) {
	foreach my $region (@{$record->{'regions'}}) {
		foreach my $orf (@{$region->{'orfs'}}) {
			$orf->{'anchor'} = $region->{'anchor'};
			$orf->{'seq_id'} = $record->{'seq_id'};
			push(@orfList, $orf);
		}
	}
}
foreach my $orf (@orfList) {
	foreach my $column (keys %$orf) {
		next if($columnHash{$column});
		push(@columnList, $column);
		$columnHash{$column} = 1;
	}
}
%columnHash = ();
foreach my $orf (@orfList) {
	foreach my $column (keys %$orf) {
		$columnHash{$column} = 1;
	}
}
@columnList = grep {$columnHash{$_}} @columnList;
@columnList = grep {$_ ne 'description'} @columnList;
print join("\t", @columnList), "\n";
foreach my $orf (@orfList) {
	print join("\t", map {defined($_) ? $_ : ''} @$orf{@columnList}), "\n";
}
