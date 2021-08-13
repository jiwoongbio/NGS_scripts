#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;

my ($summFile, $sampleNameFile) = @ARGV;
my $summ = '';
{
	open(my $reader, $summFile);
	$summ .= $_ while(<$reader>);
	close($reader);
}
$summ = decode_json($summ);
my %sampleNameHash = ();
{
	open(my $reader, $sampleNameFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($sample, $name) = split(/\t/, $line, -1);
		$sampleNameHash{$sample} = $name;
	}
	close($reader);
}
my $params;
foreach my $key (sort keys %{$summ->{'params'}}) {
	if(ref($summ->{'params'}->{$key}) eq 'HASH') {
		foreach my $sample (sort keys %{$summ->{'params'}->{$key}}) {
			$params->{$key}->{$sampleNameHash{$sample}} = $summ->{'params'}->{$key}->{$sample};
		}
	} elsif(ref($summ->{'params'}->{$key}) eq 'ARRAY') {
		foreach my $sample (@{$summ->{'params'}->{$key}}) {
			push(@{$params->{$key}}, $sampleNameHash{$sample});
		}
	}
}
$summ->{'params'} = $params;
print encode_json($summ);
