#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::Graphics;
use Bio::SeqFeature::Generic;

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'k=s' => \(my $keyStyle = 'between'),
	'w=i' => \(my $width = 800),
	'l=i' => \(my $pad_left = 10),
	'r=i' => \(my $pad_right = 10),
	's' => \(my $stranded = ''),
	'n=i' => \(my $nameIndex = ''),
);
my ($regionFile, $chromosome, $start, $end) = @ARGV;

my @featureList = ();
open(my $reader, $regionFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	if($tokenList[0] eq $chromosome && $tokenList[1] <= $end && $start <= $tokenList[2]) {
		my ($chromosome, $start, $end) = @tokenList[0, 1, 2];
		my $feature = Bio::SeqFeature::Generic->new(-start => $start, -end => $end);
		if($stranded ne '') {
			my $strand = $tokenList[3];
			$feature->strand( 1) if($strand eq '+');
			$feature->strand(-1) if($strand eq '-');
		}
		if($nameIndex ne '') {
			$feature->display_name($tokenList[$nameIndex]);
		}
		push(@featureList, $feature);
	}
}
close($reader);
@featureList = sort {$a->start <=> $b->start || $a->end <=> $b->end} @featureList;

my @featureListList = ();
foreach my $feature (@featureList) {
	my $index = 0;
	for(; $index < scalar(@featureListList); $index += 1) {
		last if($featureListList[$index]->[-1]->end < $feature->start);
	}
	push(@{$featureListList[$index]}, $feature);
}

my $panel = Bio::Graphics::Panel->new(-start => $start, -end => $end, -key_style => $keyStyle, -width => $width, -pad_left => $pad_left, -pad_right => $pad_right);
$panel->add_track(Bio::SeqFeature::Generic->new(-start => $start, -end => $end), -glyph => 'arrow', -tick => 2);
foreach my $featureList (@featureListList) {
	$panel->add_track(generic => $featureList, -glyph => 'generic', -stranded => 1, -label => 1);
}
print $panel->png;
exit 0;
