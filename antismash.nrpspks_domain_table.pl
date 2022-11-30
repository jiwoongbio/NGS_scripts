#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;

my ($antismashOutputDirectory) = @ARGV;
my @columnList = ('anchor', 'seq_id', 'id', 'type', 'start', 'end', 'predictions', 'napdoslink', 'blastlink', 'sequence', 'dna_sequence', 'abbreviation', 'html_class');
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
my %anchorSeqIdHash = ();
foreach my $record (@{$json->{'recordData'}}) {
	foreach my $region (@{$record->{'regions'}}) {
		$anchorSeqIdHash{$region->{'anchor'}} = $record->{'seq_id'};
	}
}
my @domainList = ();
foreach my $anchor (sort keys %{$json->{'details_data'}->{'nrpspks'}}) {
	foreach my $orf (@{$json->{'details_data'}->{'nrpspks'}->{$anchor}->{'orfs'}}) {
		foreach my $domain (@{$orf->{'domains'}}) {
			$domain->{'anchor'} = $anchor;
			$domain->{'seq_id'} = $anchorSeqIdHash{$anchor};
			$domain->{'id'} = $orf->{'id'};
			$domain->{'predictions'} = join(' ', map {@$_} @{$domain->{'predictions'}});
			push(@domainList, $domain);
		}
	}
}

foreach my $domain (@domainList) {
	foreach my $column (keys %$domain) {
		next if($columnHash{$column});
		push(@columnList, $column);
		$columnHash{$column} = 1;
	}
}
%columnHash = ();
foreach my $domain (@domainList) {
	foreach my $column (keys %$domain) {
		$columnHash{$column} = 1;
	}
}
@columnList = grep {$columnHash{$_}} @columnList;
@columnList = grep {$_ ne 'blastlink' && $_ ne 'sequence' && $_ ne 'dna_sequence' && $_ ne 'napdoslink'} @columnList;
print join("\t", @columnList), "\n";
foreach my $domain (@domainList) {
	print join("\t", map {defined($_) ? $_ : ''} @$domain{@columnList}), "\n";
}
