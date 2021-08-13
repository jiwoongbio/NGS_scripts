#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use JSON;
use Statistics::R;
use List::Util qw(sum max);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'n' => \(my $addNormal = ''),
	'C' => \(my $useCancerCellFraction = ''),
	'w=i' => \(my $width = 8),
	'h=i' => \(my $height = 8),
	'c=s' => \(my $colors = '#880000,#008800,#000088,#888800,#880088,#008888'),
	's=i' => \(my $maximumVertexSize = 32),
);
my ($summFile, $treeIndex, $plotFile) = @ARGV;
my $summ = '';
{
	open(my $reader, $summFile);
	$summ .= $_ while(<$reader>);
	close($reader);
}
$summ = decode_json($summ);
my @edgeList = ();
foreach my $vertex1 (sort {$a <=> $b} keys %{$summ->{'trees'}->{$treeIndex}->{'structure'}}) {
	foreach my $vertex2 (sort {$a <=> $b} @{$summ->{'trees'}->{$treeIndex}->{'structure'}->{$vertex1}}) {
		push(@edgeList, [sort {$a <=> $b} ($vertex1, $vertex2)]);
	}
}
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
my %vertexSizeHash = ();
foreach my $population (@populationList) {
	$vertexSizeHash{$population} = sum(0, map {defined($_) ? $_ : 0} @{$populationSampleFractionHash{$population}}{@sampleList});
}
#foreach my $population (@populationList) {
#	$vertexSizeHash{$population} -= sum(0, map {$vertexSizeHash{$_->[1]}} grep {$_->[0] == $population} @edgeList);
#}
my @vertexSizeList = @vertexSizeHash{@populationList};
@vertexSizeList = map {$_ / max(@vertexSizeList) * $maximumVertexSize} @vertexSizeList;
{
	my $R = Statistics::R->new();
	$R->run('edges <- c()');
	foreach my $index (0 .. $#edgeList) {
		$R->set(sprintf('edges[%d]', $index * 2 + 1), $edgeList[$index]->[0]);
		$R->set(sprintf('edges[%d]', $index * 2 + 2), $edgeList[$index]->[1]);
	}
	$R->run('edges <- factor(edges)');
	$R->run('names <- levels(edges)');
	$R->run('edges <- as.numeric(edges)');

	$R->run('library(igraph)');
	$R->run('set.seed(0)');
	$R->run('g <- graph(edges)');

	$R->run('vertex.size <- c()');
	foreach my $index (0 .. $#vertexSizeList) {
		$R->set(sprintf('vertex.size[%d]', $index + 1), $vertexSizeList[$index]);
	}
	$R->run(sprintf('pdf(file = "%s", width = %d, height = %d)', $plotFile, $width, $height));
	if(scalar(@sampleList) > 1) {
		$R->run('vertex.pie <- list()');
		foreach my $population (@populationList) {
			$R->run(sprintf('vertex.pie[["%s"]] <- c(%s)', $population, join(',', map {defined($_) ? $_ : 0} @{$populationSampleFractionHash{$population}}{@sampleList})));
		}
		$R->run(sprintf('color <- c(%s)', join(',', map {"\"$_\""} split(/,/, $colors))));
		$R->run('plot.igraph(g, vertex.shape = "pie", vertex.pie = vertex.pie, vertex.pie.color = list(color), vertex.size = vertex.size, vertex.label = names, vertex.label.dist = (vertex.size ^ 0.5) / 2, vertex.label.degree = pi, layout = layout.reingold.tilford(g, root = 1))');
		$R->run(sprintf('sample <- c(%s)', join(',', map {"\"$_\""} @sampleList)));
		$R->run('legend("topleft", legend = sample, col = color, pch = 19)');
	} else {
		$R->run('plot.igraph(g, vertex.size = vertex.size, vertex.label = names, vertex.label.dist = (vertex.size ^ 0.5) / 2, vertex.label.degree = pi, layout = layout.reingold.tilford(g, root = 1))');
	}
	$R->run('dev.off()');
	$R->stop();
}
