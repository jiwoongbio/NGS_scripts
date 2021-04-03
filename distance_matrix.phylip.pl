#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($distanceMatrixFile, $sampleNameFile) = @ARGV;
my %sampleNameHash = ();
if(defined($sampleNameFile)) {
	open(my $reader, $sampleNameFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($sample, $name) = split(/\t/, $line, -1);
		$sampleNameHash{$sample} = $name;
	}
	close($reader);
}
open(my $reader, $distanceMatrixFile);
chomp(my $line = <$reader>);
(undef, my @sampleList) = split(/\t/, $line, -1);
my %sampleSampleDistanceHash = ();
while(my $line = <$reader>) {
	chomp($line);
	my %sampleDistanceHash = ();
	(my $sample, @sampleDistanceHash{@sampleList}) = split(/\t/, $line, -1);
	$sampleSampleDistanceHash{$sample} = \%sampleDistanceHash;
}
close($reader);
foreach my $sample (@sampleList) {
	$sampleSampleDistanceHash{$sample}->{$sample} = 0;
}
print scalar(@sampleList), "\n";
foreach my $sample (@sampleList) {
	my $name = $sampleNameHash{$sample};
	$name = $sample unless(defined($name));
	my $sampleDistanceHash = $sampleSampleDistanceHash{$sample};
	print join(' ', sprintf('%-10s', $name), @$sampleDistanceHash{@sampleList}), "\n";
}
