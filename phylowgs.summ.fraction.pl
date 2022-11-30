#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'n' => \(my $addNormal = ''),
	'c' => \(my $useCancerCellFraction = ''),
);
my ($summFile, $treeIndex) = @ARGV;
my $summ = '';
{
	open(my $reader, $summFile);
	$summ .= $_ while(<$reader>);
	close($reader);
}
$summ = decode_json($summ);
my @sampleList = @{$summ->{'params'}->{'samples'}};
my @populationList = sort {$a <=> $b} keys %{$summ->{'trees'}->{$treeIndex}->{'populations'}};
my %populationSampleFractionHash = ();
foreach my $population (@populationList) {
	my @fractionList = @{$summ->{'trees'}->{$treeIndex}->{'populations'}->{$population}->{'cellular_prevalence'}};
	foreach my $index (0 .. $#sampleList) {
		$populationSampleFractionHash{$population}->{$sampleList[$index]} = $fractionList[$index];
	}
}
if($useCancerCellFraction) {
	foreach my $sample (@sampleList) {
		my $cellularPrevalence = $populationSampleFractionHash{1}->{$sample};
		$populationSampleFractionHash{$_}->{$sample} /= $cellularPrevalence foreach(grep {$_ > 0} @populationList);
	}
}
if($addNormal) {
	print join("\t", 'population', 'normal', @sampleList), "\n";
	foreach my $population (@populationList) {
		if($population > 0) {
			print join("\t", $population, 0, @{$populationSampleFractionHash{$population}}{@sampleList}), "\n";
		} else {
			print join("\t", $population, scalar(@sampleList), (0) x scalar(@sampleList)), "\n";
		}
	}
} else {
	print join("\t", 'population', @sampleList), "\n";
	foreach my $population (@populationList) {
		print join("\t", $population, @{$populationSampleFractionHash{$population}}{@sampleList}), "\n";
	}
}
