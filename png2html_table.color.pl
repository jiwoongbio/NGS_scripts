#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use GD;

my ($pngFile) = @ARGV;
my $image = GD::Image->newFromPng($pngFile);

my ($width, $height) = $image->getBounds;

my %colorCountHash = ();
foreach my $x (0 .. $width) {
	foreach my $y (0 .. $height) {
		my ($index) = $image->getPixel($x, $y);
		my ($r, $g, $b) = $image->rgb($index);
		(my $color = sprintf('#%02x%02x%02x', int($r), int($g), int($b))) =~ tr/a-z/A-Z/;
		$colorCountHash{$color} += 1;
	}
}
my @colorList = sort {$colorCountHash{$b} <=> $colorCountHash{$a} || $a cmp $b} keys %colorCountHash;

print "<table>\n";
foreach my $color (@colorList) {
	print "<tr><td style=\"background-color: $color;\">&nbsp;</td><td>$color</td><td>$colorCountHash{$color}</td></tr>\n";
}
print "</table>\n";
