#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum max min);
use GD;
use Getopt::Long qw(:config no_ignore_case);

my @colorList = ();
my @ignoreColorList = ();
GetOptions(
	'n=i' => \(my $numberOfColumns = 0),
	'r=f' => \(my $widthHeightRatio = 1),
	'c=s' => \@colorList,
	'C=s' => \@ignoreColorList,
	'M' => \(my $useMajorColor),
);
my @rgbList = ();
foreach(@colorList) {
	my $color = $_;
	$color =~ s/^#//;
	my $r = hex(substr($color, 0, 2, ''));
	my $g = hex(substr($color, 0, 2, ''));
	my $b = hex(substr($color, 0, 2, ''));
	die if($color ne '');
	push(@rgbList, [$r, $g, $b]);
}
my %ignoreColorHash = ();
$ignoreColorHash{$_} = 1 foreach(@ignoreColorList);

my ($pngFile) = @ARGV;
my $image = GD::Image->newFromPng($pngFile);

my ($width, $height) = $image->getBounds;
my $denominatorWidth = ($numberOfColumns == 0) ? 1 : ($width + 1) / $numberOfColumns;
my $denominatorHeight = $denominatorWidth / $widthHeightRatio;

my $rgb;
foreach my $x (0 .. $width) {
	foreach my $y (0 .. $height) {
		my ($index) = $image->getPixel($x, $y);
		my ($r, $g, $b) = $image->rgb($index);
		my $col = int($x / $denominatorWidth);
		my $row = int($y / $denominatorHeight);
		push(@{$rgb->[$row]->[$col]->{'r'}}, $r);
		push(@{$rgb->[$row]->[$col]->{'g'}}, $g);
		push(@{$rgb->[$row]->[$col]->{'b'}}, $b);
		(my $color = sprintf('#%02x%02x%02x', int($r), int($g), int($b))) =~ tr/a-z/A-Z/;
	}
}

print "<table>\n";
foreach my $row (0 .. scalar(@$rgb) - 1) {
	my @tdList = ();
	foreach my $col (0 .. scalar(@{$rgb->[$row]}) - 1) {
		my @rList = @{$rgb->[$row]->[$col]->{'r'}};
		my @gList = @{$rgb->[$row]->[$col]->{'g'}};
		my @bList = @{$rgb->[$row]->[$col]->{'b'}};
		if($useMajorColor) {
			my %rgbCountHash = ();
			foreach my $index (0 .. max($#rList, $#gList, $#bList)) {
				$rgbCountHash{join("\t", $rList[$index], $gList[$index], $bList[$index])} += 1;
			}
			@rList = ();
			@gList = ();
			@bList = ();
			my $count = max(values %rgbCountHash);
			foreach my $rgb (keys %rgbCountHash) {
				if($rgbCountHash{$rgb} == $count) {
					my ($r, $g, $b) = split(/\t/, $rgb, -1);
					push(@rList, $r);
					push(@gList, $g);
					push(@bList, $b);
				}
			}
		}
		my $r = mean(@rList);
		my $g = mean(@gList);
		my $b = mean(@bList);
		my $color = getColor($r, $g, $b);
		if($ignoreColorHash{$color}) {
			push(@tdList, "<td></td>");
		} else {
			push(@tdList, "<td style=\"background-color: $color;\"></td>");
		}
	}
	print join('', '<tr>', @tdList, '</tr>'), "\n";
}
print "</table>\n";

sub getColor {
	my ($r, $g, $b) = @_;
	if(@rgbList) {
		my @rgbDistanceList = map {[@$_, (($_->[0] - $r) ** 2 + ($_->[1] - $g) ** 2 + ($_->[2] - $b) ** 2) ** 0.5]} @rgbList;
		my $distance = min(map {$_->[3]} @rgbDistanceList);
		@rgbDistanceList = grep {$_->[3] == $distance} @rgbDistanceList;
		($r, $g, $b) = @{$rgbDistanceList[0]};
	}
	(my $color = sprintf('#%02x%02x%02x', int($r), int($g), int($b))) =~ tr/a-z/A-Z/;
	return $color;
}

sub mean {
	my @valueList = @_;
	return sum(@valueList) / scalar(@valueList);
}
